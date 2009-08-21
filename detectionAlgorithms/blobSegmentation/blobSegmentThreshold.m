function maskBlobs = blobSegmentThreshold(image,minSize,plotRes)
%BLOBSEGMENTTHRESHOLD segments blobs in 2D images via combination of Otsu and Rosin thresholding
%
%SYNOPSIS maskBlobs = blobSegmentThreshold(image,minSize,plotRes)
%
%INPUT  image     : 2D image to be segmented.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%
%OUTPUT maskBlobs : Mask of blobs. 1 inside blobs, 0 outside.
%
%REMARKS While the code is in principle general, it has been extensively tested only on focal adhesions.
%
%Khuloud Jaqaman August 2009

%% Output
maskBlobs = [];

%% Input

%check number of input arguments
if nargin < 1
    disp('Please enter at least image to be segmented');
    return
end

%minimum blob size
if nargin < 2 || isempty(minSize)
    minSize = 20;
end

%plot results
if nargin < 3 || isempty(plotRes)
    plotRes = 0;
end

%% Segmentation

%make sure that image is in double format
image = double(image);

%remove noise by filtering image with a Gaussian whose sigma = 1 pixel
imageFiltered = Gauss2D(image,1,1);

%estimate background by filtering image with a Gaussian whose sigma = 10
%pixels
imageBackground = Gauss2D(image,10,1);

%calculate noise-filtered and background-subtracted image
imageFilteredMinusBackground = imageFiltered - imageBackground;

%enhance features by performing a maximum filter
[sizeX,sizeY] = size(imageFilteredMinusBackground);
imageDilated = imageFilteredMinusBackground;
imageTmp(:,:,1) = imageDilated;
imageTmp(:,:,2) = [zeros(1,sizeY); imageDilated(2:end,:)];
imageTmp(:,:,3) = [imageDilated(1:end-1,:); zeros(1,sizeY)];
imageTmp(:,:,4) = [zeros(sizeX,1) imageDilated(:,2:end)];
imageTmp(:,:,5) = [imageDilated(:,1:end-1) zeros(sizeX,1)];
imageTmp(:,:,6) = [zeros(1,sizeY); [zeros(sizeX-1,1) imageDilated(2:end,2:end)]];
imageTmp(:,:,7) = [zeros(1,sizeY); [imageDilated(2:end,1:end-1) zeros(sizeX-1,1)]];
imageTmp(:,:,8) = [[zeros(sizeX-1,1) imageDilated(1:end-1,2:end)]; zeros(1,sizeY)];
imageTmp(:,:,9) = [[imageDilated(1:end-1,1:end-1) zeros(sizeX-1,1)]; zeros(1,sizeY)];
imageDilated = max(imageTmp,[],3);

%get minumum and maximum pixel values in image
minSignal = min(imageDilated(:));
maxSignal = max(imageDilated(:));

%normalize image between 0 and 1
imageDilatedNorm = (imageDilated - minSignal)/(maxSignal - minSignal);

%estimate the intensity level to use for thresholding the image
level1 = graythresh(imageDilatedNorm); %Otsu
[dummy, level2] = cutFirstHistMode(imageDilatedNorm(:),0); %Rosin
level = 0.33333*level1 + 0.66667*level2;

%threshold the image
imageThresholded = im2bw(imageDilatedNorm,level);
imageThresholded = double(imageThresholded);

%fill holes in thresholded image to make continuous blobs
imageThresholdedFilled = imfill(imageThresholded,'holes');

%find pixels belonging to blobs
[indxBlob,numBlob] = bwlabel(imageThresholdedFilled);

%find average background intensity
% intensityBG = mean(imageDilatedNorm(indxFA==0));

%go over blobs and find their size and average intensity
%discard those in which there is little confidence (generally, small size
%and low intensity) - although intensity criterion is not looked at any more
[sizeBlob,intensityBlob] = deal(zeros(numBlob,1));
for iBlob = 1 : numBlob

    %find pixel indices belonging to this blob
    indxBlobCurrent = find(indxBlob == iBlob);

    %find number of pixels belonging to this Blob
    sizeBlobCurrent = length(indxBlobCurrent);

    %get average intensity of blob
    intensityBlobCurrent = mean(imageDilatedNorm(indxBlobCurrent));

    %store size and intensity
    sizeBlob(iBlob) = sizeBlobCurrent;
    intensityBlob(iBlob) = intensityBlobCurrent;

    %discard if blob size is smaller than the minimum allowed
    if sizeBlobCurrent < minSize
        indxBlob(indxBlobCurrent) = 0;
    end 

end

%output final blob mask
maskBlobs = double(indxBlob > 0);

%% Plotting

if plotRes

    %get the blob edges from the final blob mask
    edgesBlobs = double(edge(maskBlobs));

    %scale the original image to be between 0 and 1
    imageScaled = (image - min(image(:))) / (max(image(:)) - min(image(:)));

    %give the edge pixels a value of zero in the original image
    imageScaled(edgesBlobs==1) = 0;

    %construct a 3-layered image to show blob edges on top of
    %original image
    image3Color = repmat(imageScaled,[1 1 3]);
    image3Color(:,:,1) = image3Color(:,:,1) + edgesBlobs;
    
    %plot image
    figure, hold on
    subplot(1,2,1)
    imshow(image,[])
    subplot(1,2,2)
    imshow(image3Color,[]);

end

%% extra stuff

%otherwise, if smaller than the intermediate size ...
%     elseif sizeFACurrent < intSize
%
%         %calculate actual intensity factor (ratio of FA intensity to
%         %background)
%         intensityFactor = intensityFACurrent / intensityBG;
%
%         %discard if intensity factor is smaller than minimum allowed for
%         %this size
%         if intensityFactor < minSizeMinIntensityFactor
%             indxFA(indxFACurrent) = 0;
%         end
%
%         %otherwise, if smaller than the maximum size with extra intensity
%         %penalty
%     elseif sizeFACurrent < maxSizeIntensityPenalty
%
%         %calculate minimum allowed intensity factor
%         minIntensityFactor = yIntersept + lineSlope * sizeFACurrent;
%
%         %calculate actual intensity factor (ratio of FA intensity to
%         %background)
%         intensityFactor = intensityFACurrent / intensityBG;
%
%         %discard if intensity factor is smaller than minimum allowed for
%         %this size
%         if intensityFactor < minIntensityFactor
%             indxFA(indxFACurrent) = 0;
%         end
%
%         %otherise ...
%     else
%
%         %calculate actual intensity factor (ratio of FA intensity to
%         %background)
%         intensityFactor = intensityFACurrent / intensityBG;
%
%         %discard if intensity factor is smaller than minimum allowed
%         %overall
%         if intensityFactor < minIntensityFactorAll
%             indxFA(indxFACurrent) = 0;
%         end



%extract focal adhesion selection criteria
% minSize = selectParam.minSize;
% minSizeMinIntensityFactor = selectParam.minSizeMinIntensityFactor;
% maxSizeIntensityPenalty = selectParam.maxSizeIntensityPenalty;
% minIntensityFactorAll = selectParam.minIntensityFactorAll;
%
% %calculate straight line equation for focal adhesion selection
% intSize = minSize + (maxSizeIntensityPenalty - minSize) / 3;
% % intSize = minSize;
% x1 = intSize;
% y1 = minSizeMinIntensityFactor;
% x2 = maxSizeIntensityPenalty;
% y2 = minIntensityFactorAll;
% lineSlope = (y1 - y2) / (x1 - x2);
% yIntersept = y1 - lineSlope * x1;

