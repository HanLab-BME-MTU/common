function getSubsequenceOfFiles(stepSize,outFName)
%getSubsequenceOfFiles: Get a subsequence of files every 'stepSize' from indexed
%                       files in 'inDir' and output it to 'outDir'.
%
% SYNOPSIS: getSubsequenceOfFiles(inDir,outDir,inFName,outFName,stepSize)
%
% INPUT:
%    inDir   : The path name of the input directory.
%    outDir  : The path name of the output directory.
%    inFName : The common prefix name of the input files from which the subsequence
%              is drawn.
%    outFName: The common prefix name for the output files. Pass [] to use the
%              same name as 'inFName'.
%    stepSize: Get files every 'stepSize'.
%
% Version: MATLAB Version 7.0.4.352 (R14) Service Pack 2.
% OS     : Linux.
%
% Author: Lin Ji, Dec. 20, 2005.

[firstFileName inDir filterIndex] = uigetfile('*.*','Pick the first file');
if isequal(firstFileName,0) || isequal(inDir,0)
   disp('User pressed cancel. No file is selected.');
   return;
end

[path,inFName,no,ext] = getFilenameBody(firstFileName);
firstFileIndex = str2num(no);

if isempty(outFName)
   outFName = inFName;
end
[filelist,index] = getNamedFiles(inDir,inFName);
selIndex = find(index>=firstFileIndex);
index    = index(selIndex);
filelist = filelist(selIndex);

numTotalFiles  = length(index);
subSequence    = [1:stepSize:numTotalFiles];
numSubseqFiles = length(subSequence);
indexForm      = sprintf('%%.%dd',length(num2str(numSubseqFiles)));

parentDir = [inDir filesep '..'];
outDir = uigetdir(parentDir,'Select output directory');

if isequal(outDir,0)
   disp('User pressed cancel. No output directory is selected.');
   return;
end

for k = 1:numSubseqFiles
   inFileName  = [inDir filesep filelist{subSequence(k)}];
   outFileName = [outDir filesep outFName sprintf(indexForm,k) ext];
   [success,msgID] = copyfile(inFileName,outFileName);
end
