function Im = subResSegment2DImageModel(segmentParams, imSize)
% Compute an image model which is the sum of sub-resolution 2D segments
% defined by params (see subResSegment2D.m for more details).
%
% Im = subResSegment2DImageModel(params, imSize)
%
% INPUT:
% segmentParams   : nx6 matrix where n is the number of segments and their
%                   parameters (see subResSegment2D() for details)
%
% imSize          : image size
%
% OUTPUT:
% Im              : the image model of size imSize.
%
% Sylvain Berlemont, 2010

[n,p] = size(segmentParams);

if p ~= 6
    error('Invalid number of segment parameters.');
end

% Generate the image segments
Im = zeros(imSize);

for i = 1:n
    xC = segmentParams(i,1);
    yC = segmentParams(i,2);
    amp = segmentParams(i,3);
    sigma = segmentParams(i,4);
    l = segmentParams(i,5);
    theta = segmentParams(i,6);

    [xRange,yRange,nzIdx] = subResSegment2DSupport(xC,yC,sigma,l,theta,imSize);

    S = subResSegment2D(xRange-xC,yRange-yC,amp,sigma,l,theta,nzIdx);
    
    Im(yRange,xRange) = Im(yRange,xRange) + S;
end