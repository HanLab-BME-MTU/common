function []=batchNucleiiProteinColoc
% batchNecleiiProteinColoc interactively collects project folder and files
% of interest to perform nucleus segmentation and protein detection by
% thresholding and generate table featuring detection rate (detected
% nucleii / all nucleii)

% Get the project folder path
pathProject = uigetdir(pwd,'Choose project folder');
[fileImgDAPI, pathImgDAPI] = uigetfile('*.tif','Select DAPI channel image of the first microscopy image set. Use Ctrl for multiselection',...
    'MultiSelect','on');
if ~iscell(fileImgDAPI)
    fileImgDAPI={fileImgDAPI};
end
[~,nameDAPI1]=fileparts(fileImgDAPI{1});
nameSlide1= strrep(nameDAPI1, 'c1', '');
[~,nameDAPIend]=fileparts(fileImgDAPI{end});
nameSlideEnd= strrep(nameDAPIend, 'c1', '');
nameSlideEnd = nameSlideEnd(end-2:end);

outputPathFiles = fullfile(pathImgDAPI, [nameSlide1 '-' nameSlideEnd]);
if ~exist(outputPathFiles,'dir')
    mkdir(outputPathFiles)
end
pathResults1 = [outputPathFiles filesep 'detectionSummary.csv'];
pathResults2 = [outputPathFiles filesep 'CD45intensity.csv'];
pathResults3 = [outputPathFiles filesep 'allData.mat'];
pathResults4 = [outputPathFiles filesep 'CD45intensityMedian.csv'];
numFiles = numel(fileImgDAPI);
nameSlide = cell(numFiles,1);
numDetectedNuc = zeros(numFiles,1);
numAllNuc = zeros(numFiles,1);
ratioDetected = zeros(numFiles,1);
meanIntenCD45 = cell(numFiles,1);
meanIntenCD45ND = cell(numFiles,1);
medianIntenCD45 = cell(numFiles,1);
medianIntenCD45ND = cell(numFiles,1);
readCD45 = input('Does CD45 Channel exist (1/0)?');
thresSignal = input('Threshold value for signal detection (default:3000 for Tra, 2000 for Bstrong)?');
if isempty(thresSignal)
    if readCD45
        thresSignal = 3000;%quantile(imgSignalControl(:),0.999999); % this gives around 2500
    else
        thresSignal = 2000;%quantile(imgSignalControl(:),0.999999); % this gives around 2500
    end
end    
%% Run nucleusProteinColoc
for ii=1:numFiles
    curPathImgDAPI = fullfile(pathImgDAPI,fileImgDAPI{ii});
    [nameSlide(ii),numDetectedNuc(ii),numAllNuc(ii),ratioDetected(ii),meanIntenCD45{ii},meanIntenCD45ND{ii},...
        medianIntenCD45{ii},medianIntenCD45ND{ii}]=...
        nucleusProteinColoc(pathProject,curPathImgDAPI,thresSignal,readCD45);
end
%% Store output in excel format
% 1. Table of yield:
A= table(nameSlide,numDetectedNuc,numAllNuc,ratioDetected);
writetable(A,pathResults1)
save(pathResults3)
B= table(nameSlide,meanIntenCD45, meanIntenCD45ND);
writetable(B,pathResults2)
C= table(nameSlide,medianIntenCD45, medianIntenCD45ND);
writetable(C,pathResults4)

