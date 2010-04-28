function [aVd,aStd,nPts,distVals] = intensityVsDistFromEdge(image,mask,distVals)
%INTENSITYVSDISTFROMEDGE analyzes how the image intensity varies with distance from the mask edge 
% 
% aVd = intensityVsDistFromEdge(image,mask)
% [aVd,aStd,nPts,distVals] = intensityVsDistFromEdge(image,mask)
% [aVd,aStd,nPts,distVals] = intensityVsDistFromEdge(image,mask,distVals)
% 
% This function analyzes how the intensity in the input image varies with
% distance from the edge of the mask. This is accomplished via the distance
% transform.
% 
% Input:
% 
%   image - The image to analyze intensity in.
% 
%   mask - The mask containing the area to analyze.
% 
%   distVals - The distance values to average between.
%              Optional. If not input, every integer distance will be
%              averaged.
% 
% 
% 
% Output: 
% 
%   aVd - A 1xM vector containing the average intensity at every distance
%   interval specified by distVals.
% 
%   aStd - A 1xM vector containing the standard deviation of intensity at every distance
%   interval specified by distVals.
%
%   nPts - A 1xM vector containing the number of pixels at each distance
%   interval specified by distVals.
%
%   distVals - The distance intervals used for averaging.
% 
% 
% 
% 
% Hunter Elliott
% 4/2010
%

%% ----------- Input ---------- %%

image = double(image);
mask = mask >0;

distX = bwdist(~mask);

if nargin < 3 || isempty(distVals)
    distVals = 1:max(distX(:));    
end

%% ----- Analysis ----- %%

aVd = arrayfun(@(x)(mean(image(mask(:) & distX(:) >= distVals(x) & distX(:) < distVals(x+1)))),1:(length(distVals)-1));
if nargout > 1
    aStd = arrayfun(@(x)(std(image(mask(:) & distX(:) >= distVals(x) & distX(:) < distVals(x+1)))),1:(length(distVals)-1));
end
if nargout > 2
    nPts = arrayfun(@(x)(numel(image(mask(:) & distX(:) >= distVals(x) & distX(:) < distVals(x+1)))),1:(length(distVals)-1));
end


