% *************************************************************************
% Kunal Jathal
% MusixMatch
% 
% FEATURE EXTRACTOR w/built-in Beat Detection & Chorus Location
%
% Name:     FeatureExtractor
%
% Description:
%
% This function implements extracts features from an audio snippet. It does
% this by basically first extracting the tempo of the snippet, and then
% superimposing a beat-grid onto it. The beat-grid is then used to window
% the snippet, and the energies of the windows are computed. In this
% fashion, the energy contour of the signal is obtained. The first order
% differential of the energy curve is then computed and used to extract the
% rate of energy flux, through an FFT calculation on the differential
% itself. The lowest FFT coefficient is returned as a feature.
% 
% Usage
% 
% FeatureExtractor takes in the fileName of the audio snippet to be
% analyzed, and a flag which indicates if the chorus should attempt to be
% located or not. Chorus location is roughly calculated through onset
% detection of the first order differential energy curve.
% *************************************************************************
function f = FeatureExtractor(fileName, playChorus)

% Load audio snippet.
[audioSnippet, sr] = wavread(fileName);

% Convert to mono & normalize audio for consistency
audioSnippet = mean(audioSnippet, 2);
audioSnippet = 0.99 * audioSnippet/max(abs(audioSnippet));

% We want to window the audio snippet by beat. Because most musical changes
% don't happen in time frames of less than 2 bars, slow the original beat
% down for windowing purposes.

beatSlowingFactor = 8; % Beat slowing factor
beatOriginalTimes = beat(audioSnippet,sr); % Original tempo
beatSlowLength = ceil(length(beatOriginalTimes)/beatSlowingFactor); % Slow Beat Length
beatSlow = zeros(1, beatSlowLength); % Slow Beat Times
countB1 = 1; % Fast Beat Counter
countB2 = 1; % Slow Beat Counter

% Slow the beat down.
while (countB1 <= length(beatOriginalTimes))
    beatSlow(countB2) = beatOriginalTimes(countB1);
    countB1 = countB1 + beatSlowingFactor;
    countB2 = countB2 + 1;
end

% Set Differential Gradient (DG) variables 
dglLoc = 1; % DG Lookahead for Chorus Location Detection
dglFFT = 1; % DG Lookahead for FFT
stftl = 128; % DG Short Term FT Length

% Initialize Energy Matrix
MAT = zeros(1, length(beatSlow));

% Extract features from windows
for bSlowCounter=1:length(beatSlow)        

    beatInterval = 0;

    % If we are at the last beat, the window extends to the end of the song
    if (bSlowCounter == length(beatSlow))
        beatInterval = (length(audioSnippet)/sr) - beatSlow(bSlowCounter);
    else
        beatInterval = beatSlow(bSlowCounter+1) - beatSlow(bSlowCounter); % in seconds
    end

    beatIntervalSamples = beatInterval * sr; % in samples

    % Window the signal appropriately
    frameStart = ceil(beatSlow(bSlowCounter) * sr);
    frameEnd = frameStart + ceil(beatIntervalSamples - 1);
    windowedSignal = audioSnippet(frameStart:frameEnd);

    % Get the RMS Energy of the window
    RMSEnergy = sqrt((1/length(windowedSignal)) * sum(windowedSignal.^2));
    MAT(bSlowCounter) = RMSEnergy;
end

% We now have a vector of the energies for all the windows. Because
% energies can fluctuate quite a bit between bars, smooth the energy
% countour out a bit by applying a simple moving average filter.
smoothingFactor = 2;
filterDenominator = (1/smoothingFactor)*ones(1,smoothingFactor);
filterNumerator = 1;

% Smooth out energy contour
MATsmooth = filter(filterDenominator,filterNumerator, MAT);

% Normalize the energy vector
MATsmoothNormalized = (MATsmooth/max(MATsmooth));

% The goal here is to use the rate of fluctuation of energy to eventually
% extract the feature needed. Calculate the first order differential of the
% energy curve to get some meaningful info about it's rate. We use two
% separate first order differential curves; one for feature extraction and
% another for potential chorus location detection.
firstOrderDifferentialFFT = zeros(1, length(MATsmoothNormalized)-dglFFT);
firstOrderDifferentialLoc = zeros(1, length(MATsmoothNormalized)-dglLoc);

for windowCount=1:length(MATsmoothNormalized)-dglFFT
    vecOne = MATsmoothNormalized(windowCount);
    vecTwo = MATsmoothNormalized(windowCount + dglFFT);    
    firstOrderDifferentialFFT(windowCount) = vecTwo - vecOne;
end

for windowCount=1:length(MATsmoothNormalized)-dglLoc
    vecOne = MATsmoothNormalized(windowCount);
    vecTwo = MATsmoothNormalized(windowCount + dglLoc);    
    firstOrderDifferentialLoc(windowCount) = vecTwo - vecOne;
    
    % We weigh the first order differential curve to compensate for
    % transitions that don't correspond to choruses. This makes chorus
    % location detection a little more precise/robust.
    meanVec = (vecTwo + vecOne) / 2;
    firstOrderDifferentialLoc(windowCount) = firstOrderDifferentialLoc(windowCount)/meanVec;
end

% Audio snippets with chorus sections tend to have more energy fluctuation
% than snippets without chorus. This means that the curve of the energy
% contour will have a higher frequency. So we can take the FFT of the first
% order differential and analyze the low frequency bin. Chorus sections
% tend to have more high frequency (i.e. less low frequency) than
% Non-chorus sections. This is sneaky, but it works.
feature = abs(fft(firstOrderDifferentialFFT, stftl))/length(firstOrderDifferentialFFT);

% Analyze the lowest bin, since we don't need all the rest.
f = feature(1);

if (playChorus)
    % Now that we have the differential, let us use it to split the song into
    % segments.

    % We will achieve this by building an onset detector, that will use the
    % mean of the contour as it's threshold. The beginning of audio is the
    % first segment. When we pass the threshold, the next segment starts, until
    % we fall below the threshold. This should roughly split the song into chorus 
    % and non-chorus segments.

    segmentLocations = 1;
    currentBeat = 2;
    belowThreshold = false;
    aboveThreshold = false;
    neutral = true;

    % Set threshold
    threshold = mean(abs(firstOrderDifferentialLoc));
    
    % Onset Detection
    for currentBeat = 1:length(firstOrderDifferentialLoc)
        if (neutral)        
            if (firstOrderDifferentialLoc(currentBeat) > threshold)
                segmentLocations = [segmentLocations currentBeat+dglLoc];
                neutral = false;
                aboveThreshold = true;
                belowThreshold = false;
            elseif (firstOrderDifferentialLoc(currentBeat) < -threshold)
                segmentLocations = [segmentLocations currentBeat+dglLoc];
                neutral = false;
                aboveThreshold = false;
                belowThreshold = true;
            end        
        elseif (aboveThreshold)
            if (firstOrderDifferentialLoc(currentBeat) < threshold)
                if (firstOrderDifferentialLoc(currentBeat) < -threshold)
                    belowThreshold = true;
                    neutral = false;
                    aboveThreshold = false;
                    segmentLocations = [segmentLocations currentBeat+dglLoc];
                else
                    neutral = true;
                    aboveThreshold = false;
                    belowThreshold = false;
                end
            end
        elseif (belowThreshold)
            if (firstOrderDifferentialLoc(currentBeat) > -threshold)
                if (firstOrderDifferentialLoc(currentBeat) > threshold)
                    belowThreshold = false;
                    neutral = false;
                    aboveThreshold = true;
                    segmentLocations = [segmentLocations currentBeat+dglLoc];
                else
                    neutral = true;
                    aboveThreshold = false;
                    belowThreshold = false;
                end
            end            
        end                    
    end

    % Now that we have the different segment locations, split them up.
    numberOfSegments = length(segmentLocations);
    segment = 0;
    currentEnergy = 0;
    highestEnergy = 0;
    chorus = 0;
    
    for splitCounter=1:numberOfSegments
        startingSample = floor(beatSlow(segmentLocations(splitCounter)) * sr);
        endingSample = 0;

        % If we are at the last segment, append the ending of the song to it.
        if (splitCounter == numberOfSegments)
            endingSample = length(audioSnippet);
        else    
            endingSample = floor((beatSlow(segmentLocations(splitCounter + 1)) * sr) - 1);
        end

        segment = audioSnippet(startingSample:endingSample);
        currentEnergy = sqrt((1/length(segment)) * sum((segment.^2)));
        
        % See which segment has the highest RMS energy. It's likely to be the
        % chorus (or very near it). Because we are looking at RMS energy over
        % the entire window, this provides some resilience to sudden energy
        % transitions (egs: build-ups) that could be otherwise mistaken for
        % choruses.
        if (currentEnergy > highestEnergy)
            highestEnergy = currentEnergy;
            chorus = segment;
        end
    end
    
    % Preview 5 seconds of the chorus. Add a small fade-out for style points.
    interval = 5 * sr;
    
    % If the segment is too short, it likely isn't the chorus, so skip it.
    if (interval > length(chorus))
        disp('Chorus location uncertain, skipping preview..');
    else
        disp('Preview...');
        envelope = [ones(1, interval-sr) linspace(1, 0, sr)];    
        chorusPreview = chorus(1:interval) .* envelope';
        player = audioplayer(chorusPreview, sr);
        playblocking(player);
    end
end

end
