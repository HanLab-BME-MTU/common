function  [level,bw_out] = thresholdOtsu_local(imageIn,level_local_radius, pace, lowerbound, varargin)
% local thresholding based on thresholding level using Otsu's method
%
% [level,bw_out] = thresholdOtsu_local(imageIn,level_local_radius, pace, showPlots)
% 
% This function selects a threshold for the input fluorescence image using
% Ostu method(function thresholdOtsu()), then find the local thresholds also
% decided by Otsu method(function thresholdOtsu()) to get a local threshold
% map. After some smoothing this threshold map is used to segment the input
% image.
% 
% Input:
% 
%   imageIn:            2D input image to be thresholded.
%   level_local_radius: the radius of local patch
%   pace:               the pace to calculate local threshold, mostly to
%                       speed up the process which is computational expensive
%   lowerbound:         the percentage of global threshold the user want to
%                       set as the lowerbound of the local thresholding
%   showPlots:          If true, a plot of the histogram and an overlay of the mask
%                       on the image will be shown. 
%
% Output:
% 
%   level - The intensity value selected for global thresholding.
%   bw_out - the output segmentation result
%

% Liya Ding, 1/2012

% this function is re-organized in the function called thresholdLocal with
% the option of 'Otsu'
[level,bw_out] = thresholdLocalSeg(imageIn,'Otsu', level_local_radius, pace, lowerbound);