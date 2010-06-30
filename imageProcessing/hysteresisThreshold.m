% Performs hysteresis thresholding.

% Francois Aguet, June 30, 2010

function out = hysteresisThreshold(image, t1, t2)

im2 = image;
im2(im2<t2) = 0;
im2(im2~=0) = 1;

[labels ncomp] = bwlabel(im2, 8);

im1 = image;
im1(im1<t1) = 0;
im1(im1~=0) = 1;

labelsAboveT1 = labels.*im1;

validLabels = unique(labelsAboveT1(:));
validLabels = validLabels(2:end);

out = zeros(size(image));
for k = 1:length(validLabels)
    out(labels==validLabels(k)) = 1;
end