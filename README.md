# musixmatch
Musixmatch Chorus Detection

Usage:

1) RUN chorusClassifierFinal.m in a MATLAB environment. Minimum version required is R2012b, but R2015a recommended.

2) No arguments are needed to this function; simply run chorusClassifierFinal. 

3) Because MATLAB doesn't natively support mp3s (audioread requires the use of external binaries), the Classifier only reads wav data. Run TurnWarningsOff.m before running chorusClassifierFinal to get rid of annoying wavread warnings. 

4) If you want, you can copy any of your own audio snippets you want to use as your dataset (both testing and training) into the 'Songs' Folder. The training and testing songs can be individually listed in the code itself under their corresponding vector variables. If you use your own songs for training/testing, make sure to also update their corresponding group numbers in the group vector variables in the code. 

5) This is a chorus detector; it detects whether a chorus exists in a specific section of audio. To additionally try and locate the position of the chorus in the audio, uncomment the FeatureExtractor line of code in chorusClassifierFinal to hear small previews of every chorus correctly detected.

Shoot me an email if you have any questions!
