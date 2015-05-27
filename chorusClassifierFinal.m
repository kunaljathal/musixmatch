% *************************************************************************
% Kunal Jathal
% MusixMatch
% 
% CHORUS CLASSIFIER
%
% Name:     chorusClassifierFinal
%
% Description:
%
% This function implements a basic chorus detection classifier. It uses a
% training set of audio snippets that contain either just a verse, or a
% verse that transitions to a chorus. It analyzes and uses the frequency of
% energy fluctuations in the audio snippets to try and feed a kNN
% classifier, which eventually tries to detect if a test sample audio
% snippet contains a chorus transition or not. In addition to implementing
% the classifier, a rough chorus location system has also been implemented.
% The chorus location system tries to identify the exact position and
% duration of the chorus and plays it back.
% 
% Usage
% 
% Call this function as you would any other MATLAB function. The list of
% training snippets, groups, and test audio snippets can be specified in
% the songVectorTrain, groupVectorTrain, and testSample fields below.
% *************************************************************************
function chorusClassifierFinal

% List out the file names of the audio snippets to be used for TRAINING here
songVectorTrain = char(...
'BabyOneMoreTime_Chorus.wav',...
'BabyOneMoreTime_NoChorus.wav',...
'Umbrella_Chorus.wav',...
'Umbrella_NoChorus.wav',...
'Tubthumping_Chorus.wav',...
'Tubthumping_NoChorus.wav'...
);

% MATLAB is a bit cumbersome with string arrays, so use numerical values
% for groups. These can be easily translated later. 1 = Chorus, 0 = No
% Chorus.
groupVectorTrain = [1, 0, 1, 0, 1, 0];

% List the file name of the audio snippet you want to TEST here. (Try some
% other test samples too: CaliforniaGirls_NoChorus.wav, Happy_Chorus.wav,
% TeenageDream_Chorus.wav, Happy_NoChorus.wav, etc.)
testSample = 'CaliforniaGirls_Chorus.wav';

songList = size(songVectorTrain);
numberOfSongs = songList(1);

% Initialize empty training feature vector
featureVectorTrain = zeros(numberOfSongs, 1);

for song=1:numberOfSongs
    % First, build the training set feature and group vectors.
    fileNameTrain = strcat(sprintf('Songs/%s', songVectorTrain(song, :)));
    feature = FeatureExtractor(fileNameTrain, false);
    featureVectorTrain(song, 1) = feature;
end

% Next, build the feature to test.
fileNameTest = sprintf('Songs/%s', testSample);
featureVectorTest = FeatureExtractor(fileNameTest, false);

class = knnclassify(featureVectorTest, featureVectorTrain, groupVectorTrain, 3);

if (class == 0)
    disp(sprintf('\nAudio snippet %s does NOT contain a chorus transition.\n', testSample));
else
    disp(sprintf('\nAudio snippet %s contains a chorus transition!\n', testSample));
    
    % Attempt to locate the chorus and play it! Optional.
    % FeatureExtractor(fileNameTest, true);
end

end


