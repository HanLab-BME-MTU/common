function maskLines = lineSegmentThreshold(image,minSize,plotRes,mask)
%LINESEGMENTTHRESHOLD ...
%
%SYNOPSIS maskLines = lineSegmentThreshold(image,minSize,plotRes)
%
%INPUT  image     : 2D image to be segmented.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       mask      : Binary mask. Optional. If not provided, the whole image
%                   domain is segmented.
%
%OUTPUT maskLines : Mask of lines. 1 inside lines, 0 outside.
%
%REMARKS While the code is in principle general, it has been extensively tested only on focal adhesions.
%
%Khuloud Jaqaman August 2009

%% Output
maskLines = [];

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

%mask
if nargin < 4 || isempty(mask)
    mask = ones(size(image));
end

if ~logical(mask)
    error('Mask must be a logical image.');
end
    
%% Segmentation

%make sure that image is in double format
image = double(image);

%remove noise by filtering image with a Gaussian whose sigma = 1 pixel
% imageFiltered = filterGauss2D(image,1);
[imageFiltered,~,NMS] = steerableDetector(image,2,2);
imageFiltered = imageFiltered .* NMS;

%crop image
imageFiltered = imageFiltered .* mask;

% %enhance features by performing a maximum filter
% [sizeX,sizeY] = size(imageFiltered);
% imageDilated = imageFiltered;
% imageTmp(:,:,1) = imageDilated;
% imageTmp(:,:,2) = [zeros(1,sizeY); imageDilated(2:end,:)];
% imageTmp(:,:,3) = [imageDilated(1:end-1,:); zeros(1,sizeY)];
% imageTmp(:,:,4) = [zeros(sizeX,1) imageDilated(:,2:end)];
% imageTmp(:,:,5) = [imageDilated(:,1:end-1) zeros(sizeX,1)];
% imageTmp(:,:,6) = [zeros(1,sizeY); [zeros(sizeX-1,1) imageDilated(2:end,2:end)]];
% imageTmp(:,:,7) = [zeros(1,sizeY); [imageDilated(2:end,1:end-1) zeros(sizeX-1,1)]];
% imageTmp(:,:,8) = [[zeros(sizeX-1,1) imageDilated(1:end-1,2:end)]; zeros(1,sizeY)];
% imageTmp(:,:,9) = [[imageDilated(1:end-1,1:end-1) zeros(sizeX-1,1)]; zeros(1,sizeY)];
% imageDilated = max(imageTmp,[],3);

imageDilated = imageFiltered;

% Find non zero values (due to masking)
nzInd = find(imageDilated);

%get minumum and maximum pixel values in image
minSignal = min(imageDilated(nzInd));
maxSignal = max(imageDilated(nzInd));

%normalize non zero value between 0 and 1
imageDilatedNorm = zeros(size(imageDilated));
imageDilatedNorm(nzInd) = (imageDilated(nzInd) - minSignal) / (maxSignal - minSignal);

%estimate the intensity level to use for thresholding the image
tmp = imageDilatedNorm(nzInd);
[dummy, level2] = cutFirstHistMode(tmp,0); %Rosin
level = min(level2,prctile(tmp,80));

%threshold the image
imageThresholded = im2bw(imageDilatedNorm,level);

%fill holes in thresholded image to make continuous lines
imageThresholdedFilled = imfill(imageThresholded,'holes');

%go over lines and remove those with a size smaller that minSize
labels = bwlabel(imageThresholdedFilled);
stats = regionprops(labels, 'Area'); %#ok<MRPBW>
idx = find([stats.Area] > minSize);

%output final blob mask
maskLines = ismember(labels, idx);

%dilate masks and erode masks
SE = strel('square',5);
maskLines = imdilate(maskLines,SE);
SE = strel('square',3);
maskLines = imerode(maskLines,SE);

%% Plotting

if plotRes

    %get the blob edges from the final blob mask
    edgesLines = double(edge(maskLines));

    %scale the original image to be between 0 and 1
    imageScaled = (image - min(image(:))) / (max(image(:)) - min(image(:)));

    %give the edge pixels a value of zero in the original image
    imageScaled(edgesLines==1) = 0;

    %construct a 3-layered image to show blob edges on top of
    %original image
    image3Color = repmat(imageScaled,[1 1 3]);
    image3Color(:,:,1) = image3Color(:,:,1) + edgesLines;
    
    %plot image
    figure, hold on
    subplot(1,2,1)
    imshow(image,[])
    subplot(1,2,2)
    imshow(image3Color,[]);

end

%% ~~~ the end ~~~
