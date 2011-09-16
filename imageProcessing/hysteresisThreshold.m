%out = hysteresisThreshold(img, t1, t2) performs hysteresis thresholding.
%
% Inputs:  img : input image
%           t1 : high threshold
%           t2 : low threshold
%
% Outputs: out : binary mask

% Francois Aguet, June 30, 2010

function out = hysteresisThreshold(img, t1, t2)

im2 = img;
im2(im2<t2) = 0;
im2(im2~=0) = 1;

labels = bwlabel(im2, 8);

im1 = img;
im1(im1<t1) = 0;
im1(im1~=0) = 1;

labelsAboveT1 = labels.*im1;

validLabels = unique(labelsAboveT1(:));
validLabels = validLabels(2:end);

out = ismember(labels, validLabels);