function mask = firstMinAfterFirstMaxSeg(thismask, closureRadius, sigmaGauss, holes)
%FIRSTMINAFTERFIRSTMAXSEG Detects cell edge by thresholding the intensity histogram at the first minimum after the highest maximum
%
%SYNOPSIS mask = firstMinAfterFirstMaxSeg(thismask, closureRadius)
%
%Input:    thisimage                       target image for segmentation
%          (optional) closureRadius        performed image closing on the
%                                          segmentation mask with a disk
%                                          with radius = closureRadius
%          (optional) holes                if holes==1 the mask might
%                                          contain holes. If holes==0 all
%                                          holes will be filled up. The
%                                          default is no holes (holes=0).
%
%Output:   mask                            segmented mask
%
%Khuloud Jaqaman, March 2010

if( nargin < 2)
	closureRadius = 0;
end

if nargin < 3 || isempty(sigmaGauss)
    sigmaGauss=1;
end

if nargin < 4 || isempty(holes)
    holes=0;
end


%make sure that image is in double format
thismask = double(thismask);

%remove noise by filtering image with a Gaussian whose sigma in [pixel] = sigmaGauss 
imageFiltered = filterGauss2D(thismask,sigmaGauss);

%get minumum and maximum pixel values in image
minSignal = min(imageFiltered(:));
maxSignal = max(imageFiltered(:));

%normalize image between 0 and 1
imageFilteredNorm = (imageFiltered - minSignal)/(maxSignal - minSignal);

%estimate the intensity level to use for thresholding the image
%both of the algorithms below give the same answer
%the second is faster, thus I will use it
% level = splitModes(imageFilteredNorm(:),[0 0]);
level = thresholdFluorescenceImage(imageFilteredNorm);

%make mask
mask = im2bw(imageFilteredNorm,level);

% find the largest mask
L = bwlabel(mask);
s  = regionprops(L, 'Area');
Allarea = [s.Area];
[tmp1, tmp2] = max(Allarea);
mask = (L == tmp2);

% perform closing operation
closureBrush = strel('disk',closureRadius);
mask = imclose(cast(mask,'double'),closureBrush);
mask = cast(mask,'logical');
if holes==0
    mask = double(imfill(double(mask),'holes'));
    %mask(1,:) = 0;
    L2 = bwlabel(mask==0);
    s2  = regionprops(L2, 'Area','PixelIdxList');

    % find the background with the highest average intensity
    % assign the rest as foreground
    numberOfBlobs = size(s2, 1);
    for k=1:numberOfBlobs
        PixelIdxList = s2(k).PixelIdxList;
        meanGL(k) = mean(thismask(PixelIdxList));
    end
    [tmp1, tmp2] = min(meanGL);
    maskbg = (L2 == tmp2);
    mask = (maskbg == 0);
end