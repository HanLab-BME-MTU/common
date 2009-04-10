function mask = otsuSeg(thismask, closureRadius)
%ipexbatchDetectCells Algorithm to detect cells in image.
%   segmentedCells = ipexbatchDetectCells(I) detects cells in the cell
%   image I and returns the result in segmentedCells.
%
%   Supports batch processing demo, ipexbatch.m,
%   ("Batch Processing of Image Files Using Distributed Computing").

%   Copyright 2005 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $ $Date: 2005/06/20 03:09:37 $

thismask = nrm(thismask,16);
level = graythresh(thismask);
mask = im2bw(thismask,level);

L = bwlabel(mask);
s  = regionprops(L, 'Area');
Allarea = [s.Area];
[tmp1, tmp2] = max(Allarea);
mask = (L == tmp2);

closureBrush = strel('disk',closureRadius);
mask = imclose(cast(mask,'double'),closureBrush);
mask = cast(mask,'logical');
mask = double(imfill(double(mask),'holes'));

L2 = bwlabel(mask==0);
s2  = regionprops(L2, 'Area','PixelIdxList');

numberOfBlobs = size(s2, 1);
for k=1:numberOfBlobs
    PixelIdxList = s2(k).PixelIdxList;
    meanGL(k) = mean(thismask(PixelIdxList));
end
[tmp1, tmp2] = min(meanGL);
maskbg = (L2 == tmp2);
mask = (maskbg == 0);
