function [Im, xRange, yRange, nzIdx] = subResSegment2DImageModel(segmentParams, imSize)
% Compute an image model which is the sum of sub-resolution 2D segments
% defined by params (see subResSegment2D.m for more details).
%
% [Im, xRange, yRange, nzIdx] = subResSegment2DImageModel(params, imSize)
%
% INPUT:
% segmentParams      nx7 matrix where n is the number of segments and their
%                    parameters (see subResSegment2D() for details)
%
% imSize             image size
%
% OUTPUT:
% Im                 the image model of size imSize.
%
% (xRange, yRange)   2 cell arrays of vectors representing the
%                    2-dimensional support of the segment in the image
%                    domain 1:size(imSize,1) x 1:size(imSize,2).
%
% nzIdx              linear indices where pixel value is not zero. These
%                    indices are local and are intended to be passed to
%                    subResSegment2D() and subResSegment2DJacobian()
%                    functions.
%
% Sylvain Berlemont, 2010

[n,p] = size(segmentParams);

if p ~= 7
    error('Invalid number of segment parameters.');
end

xRange = cell(n, 1);
yRange = cell(n, 1);
nzIdx = cell(n,1);

% Generate the image segments
Im = zeros(imSize);

for i = 1:n
    xC = segmentParams(i,1);
    yC = segmentParams(i,2);
    amp = segmentParams(i,3);
    sigma = segmentParams(i,4);
    l = segmentParams(i,5);
    theta = segmentParams(i,6);
    bg = segmentParams(i,7);

    [xRange{i},yRange{i},nzIdx{i}] = subResSegment2DSupport(xC,yC,sigma,l,theta,imSize);

    S = subResSegment2D(xRange-xC,yRange-yC,amp,sigma,l,theta,bg,nzIdx{i});
    
    Im(yRange{i},xRange{i}) = Im(yRange{i},xRange{i}) + S;
end