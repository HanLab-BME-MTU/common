function []=batchNecleiiProteinColoc
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
numFiles = numel(fileImgDAPI);
nameSlide = cell(numFiles,1);
numDetectedNuc = zeros(numFiles,1);
numAllNuc = zeros(numFiles,1);
ratioDetected = zeros(numFiles,1);
meanIntenCD45 = cell(numFiles,1);
meanIntenCD45ND = cell(numFiles,1);
readCD45 = input('Read CD45 Channel (1/0)?');
%% Run nucleusProteinColoc
for ii=1:numFiles
    curPathImgDAPI = fullfile(pathImgDAPI,fileImgDAPI{ii});
    [nameSlide(ii),numDetectedNuc(ii),numAllNuc(ii),ratioDetected(ii),meanIntenCD45{ii},meanIntenCD45ND{ii}]=...
        nucleusProteinColoc(pathProject,curPathImgDAPI,readCD45);
end
%% Store output in excel format
% 1. Table of yield:
A= table(nameSlide,numDetectedNuc,numAllNuc,ratioDetected);
writetable(A,pathResults1)
save(pathResults3)
B= table(nameSlide,meanIntenCD45, meanIntenCD45ND);
writetable(B,pathResults2)

