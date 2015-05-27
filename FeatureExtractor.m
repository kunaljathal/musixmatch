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

% Load audio snippet
[audioSnippet, sr] = wavread(fileName);

% Convert to mono
audioSnippet = mean(audioSnippet, 2);

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

% Set the Differential Gradient Lookahead
dgl = 1;

% If you want to hear the beats played on top of the audio, uncomment here
% db = mkblips(bSlow,sr,length(d));
% soundsc(d+db,sr);

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

    % Get the FFT ready.
    fftSize = nextpow2(beatIntervalSamples);
    fftLength = 2^fftSize;

    frameStart = ceil(beatSlow(bSlowCounter) * sr);
    frameEnd = frameStart + ceil(beatIntervalSamples - 1);
    windowedSignal = audioSnippet(frameStart:frameEnd);

    % Compute the FFT and get the Energy of the signal
    fftWindowedSignal = abs(fft(windowedSignal, fftLength));    
    fftWindowedSignal = fftWindowedSignal.^2;
    fftWindowedSignal = fftWindowedSignal / fftLength;
    MAT(bSlowCounter) = sum(fftWindowedSignal, 1);

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
MATsmoothNormalized = (MATsmooth/max(MATsmooth)) * 1000;

% The goal here is to use the rate of fluctuation of energy to eventually
% extract the feature needed. Calculate the first order differential of the
% energy curve to get some meaningful info about it's rate.
firstOrderDifferential = zeros(1, length(MATsmoothNormalized)-dgl);

for windowCount=1:length(MATsmoothNormalized)-dgl
    vecOne = MATsmoothNormalized(windowCount);
    vecTwo = MATsmoothNormalized(windowCount + dgl);    
    firstOrderDifferential(windowCount) = vecTwo - vecOne;    
end

% Audio snippets with chorus sections tend to have more energy fluctuation
% than snippets without chorus. This means that the curve of the energy
% contour will have a higher frequency. So we can take the FFT of the first
% order differential and analyze the low frequency bin. Chorus sections
% tend to have more high frequency (i.e. less low frequency) than
% Non-chorus sections. This is sneaky, but it seems to work...
feature = abs(fft(firstOrderDifferential, 128))/length(firstOrderDifferential);

% Analyze the lowest bin, since we don't need all the rest.
f = feature(1);

if (playChorus)
    % Now that we have the differential, let us use it to split the song into
    % segments.

    % We will achieve this by building an onset detector, that will use the
    % mean of the contour as it's threshold. The beginning of audio is the
    % first segment. When we pass the threshold, the next segment starts, until
    % we fall below the threshold. This should split the song into chorus and
    % non-chorus segments.... roughly.

    threshold = mean(abs(firstOrderDifferential));

    segmentLocations = 1;
    currentBeat = 2;
    belowThreshold = false;
    aboveThreshold = false;
    neutral = true;

    for currentBeat = 1:length(firstOrderDifferential)
        if (neutral)        
            if (firstOrderDifferential(currentBeat) > threshold)
                segmentLocations = [segmentLocations currentBeat+dgl];
                neutral = false;
                aboveThreshold = true;
                belowThreshold = false;
            elseif (firstOrderDifferential(currentBeat) < -threshold)
                segmentLocations = [segmentLocations currentBeat+dgl];
                neutral = false;
                aboveThreshold = false;
                belowThreshold = true;
            end        
        elseif (aboveThreshold)
            if (firstOrderDifferential(currentBeat) < threshold)
                if (firstOrderDifferential(currentBeat) < -threshold)
                    belowThreshold = true;
                    neutral = false;
                    aboveThreshold = false;
                    segmentLocations = [segmentLocations currentBeat+dgl];
                else
                    neutral = true;
                    aboveThreshold = false;
                    belowThreshold = false;
                end
            end
        elseif (belowThreshold)
            if (firstOrderDifferential(currentBeat) > -threshold)
                if (firstOrderDifferential(currentBeat) > threshold)
                    belowThreshold = false;
                    neutral = false;
                    aboveThreshold = true;
                    segmentLocations = [segmentLocations currentBeat+dgl];
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
    index = 0;
    
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
        
        % See which segment has the highest RMS energy. Hopefully its the chorus!
        if (currentEnergy > highestEnergy)
            highestEnergy = currentEnergy;
            chorus = segment;
            index = splitCounter;
        end
    end
    
    % Play the segment with the most energy...
    soundsc(chorus, sr);
end

end
