function [Im, xRange, yRange, nzIdx] = segment2DImageModel(segmentParams, avgBkg, imSize)
% Compute an image model which is the sum of sub-resolution 2D segments
% defined by params (see segment2D.m for more details).
%
% [Im, xRange, yRange, nzIdx] = segment2DImageModel(segmentParams, imSize)
%
% INPUT:
% segmentParams      nx6 matrix where n is the number of segments and their
%                    parameters (see segment2D() for details)
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
%                    segment2D() and segment2DJacobian() functions.
%
% Sylvain Berlemont, 2010

[n,p] = size(segmentParams);

if p ~= 6
    error('Invalid number of segment parameters.');
end

xRange = cell(n, 1);
yRange = cell(n, 1);
nzIdx = cell(n,1);

% Generate the image segments
Im = ones(imSize) * avgBkg;

for i = 1:n
    x = segmentParams(i,1);
    y = segmentParams(i,2);
    amp = segmentParams(i,3);
    sigma = segmentParams(i,4);
    l = segmentParams(i,5);
    theta = segmentParams(i,6);

    [xRange{i},yRange{i},nzIdx{i}] = segment2DSupport(x,y,sigma,l,theta,imSize);

    S = segment2D(x,y,amp,sigma,l,theta,xRange{i},yRange{i},nzIdx{i});
    
    Im(yRange{i},xRange{i}) = Im(yRange{i},xRange{i}) + S;
end