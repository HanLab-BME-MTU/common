function Im = subResSegment2DImageModel(segmentParams, sigmaPSF, imSize)
% Compute an image model which is the sum of sub-resolution 2D segments
% defined by params (see subResSegment2D.m for more details).
%
% Im = subResSegment2DImageModel(params, sigmaPSF, imSize)
%
% INPUT:
% segmentParams   : nx5 matrix where n is the number of segments and their
%                 parameters, i.e. xC, yC, A, l, t are stored column-wise.
%
% sigmaPSF        : half width of the gaussian PSF model
%
% imSize          : image size
%
% OUTPUT:
% Im              : the image model of size imSize.
%
% Sylvain Berlemont, 2010

[n,p] = size(params);

if p ~= 5
    error('Invalid number of segment parameters.');
end

% Generate the image segments
Im = zeros(imSize);

for i = 1:n
    xC = segmentParams(i,1);
    yC = segmentParams(i,2);
    A = segmentParams(i,3);
    l = segmentParams(i,4);
    t = segmentParams(i,5);

    [xRange,yRange,nzIdx] = subResSegment2DSupport(xC,yC,sigmaPSF,l,t,imSize);

    S = subResSegment2D(xRange-xC,yRange-yC,A,sigmaPSF,l,t,nzIdx);
    
    Im(yRange,xRange) = Im(yRange,xRange) + S;
end