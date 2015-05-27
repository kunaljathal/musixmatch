% Turn a couple of warnings off. MATLAB whines about using audioread instead of
% wavread, but audioread requires the use of external binaries to run. So
% suppress the warning for now. Similar logic for knnclassify.
warning('off', 'MATLAB:audiovideo:wavread:functionToBeRemoved')
warning('off', 'bioinfo:knnclassify:incompatibility')