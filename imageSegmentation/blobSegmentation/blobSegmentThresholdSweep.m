function blobSegmentThresholdSweep(image,minSize,plotRes,plotName)
%BLOBSEGMENTTHRESHOLDSWEEP is a wrapper for blobSegmentThresholdGeneral
%
%SYNOPSIS blobSegmentThresholdSweep(image,minSize,plotRes)
%
%INPUT  image     : 2D image to be segmented.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       plotName  : Name for plots.
%
%OUTPUT figures from blobSegmentThresholdGeneral, saved in a directory
%       chosen by the user
%
%Khuloud Jaqaman May 2011

thresholdMethod = {'otsu','rosin','minmax'};
filterNoise = [0 1];
filterBackground = [0 1];

saveDir = uigetdir;

for k = 1 : 3
    for j = 1 : 2
        for i = 1 : 2
            blobSegmentThresholdGeneral(image,thresholdMethod{k},...
                filterNoise(j),filterBackground(i),minSize,plotRes,...
                plotName);
        end
    end
end

