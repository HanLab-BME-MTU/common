function [aVd,aStd,nPts,distVals] = intensityVsDistFromEdge(image,mask,distVals)
%INTENSITYVSDISTFROMEDGE analyzes how the image intensity varies with distance from the mask edge 
% 
% aVd = intensityVsDistFromEdge(image,mask)
% [aVd,aStd,nPts,distVals] = intensityVsDistFromEdge(image,mask)
% [aVd,aStd,nPts,distVals] = intensityVsDistFromEdge(image,mask,distVals)
% 
% This function analyzes how the intensity in the input image varies with
% distance from the edge of the mask. This is accomplished via the distance
% transform. Positive distances are inside the mask, negative distances are
% outside the mask.
% 
% Input:
% 
%   image - The image to analyze intensity in.
% 
%   mask - The mask containing the area to analyze.
% 
%   distVals - The distance values to average between.
%              Optional. If not input, every positive integer distance will
%              be averaged.
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
% Hunter Elliott
% 4/2010
%

%% ----------- Input ---------- %%

image = double(image);
mask = mask > 0;

distX = bwdist(~mask);

if nargin < 3 || isempty(distVals)
    distVals = 0:max(distX(:));    
else
    distVals = distVals(1:find(distVals<=max(distX(:)),1,'last'));
end

%Separate the positive and negative distances as these will use different
%distance transforms.
posDistVals = distVals(distVals>=0);
negDistVals = [distVals(distVals<0) min(posDistVals)];

if numel(negDistVals) > 1
    distX2 = -bwdist(mask);
    aVd = arrayfun(@(x)(mean(image(distX2(:) ~=0 & distX2(:) >= negDistVals(x) & distX2(:) < negDistVals(x+1)))),1:(length(negDistVals)-1));
else
    aVd = [];
end

%% ----- Analysis ----- %%


aVd = [aVd arrayfun(@(x)(mean(image(distX(:) ~= 0 & distX(:) > posDistVals(x) & distX(:) <= posDistVals(x+1)))),1:(length(posDistVals)-1))];

if nargout > 1
    if numel(negDistVals) > 1
        aStd = arrayfun(@(x)(std(image(distX2(:) ~=0 & distX2(:) >= negDistVals(x) & distX2(:) < negDistVals(x+1)))),1:(length(negDistVals)-1));
    else
        aStd = [];
    end
    aStd = [aStd arrayfun(@(x)(std(image(distX(:) ~= 0 & distX(:) > posDistVals(x) & distX(:) <= posDistVals(x+1)))),1:(length(posDistVals)-1))];
end
if nargout > 2
    if numel(negDistVals) > 1
        nPts = arrayfun(@(x)(numel(image(distX2(:) ~=0 & distX2(:) >= negDistVals(x) & distX2(:) < negDistVals(x+1)))),1:(length(negDistVals)-1));
    else
        nPts = [];
    end
    nPts = [nPts arrayfun(@(x)(numel(image(distX(:) ~= 0 & distX(:) > posDistVals(x) & distX(:) <= posDistVals(x+1)))),1:(length(posDistVals)-1))];
end


