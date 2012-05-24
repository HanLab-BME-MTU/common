function maskBlobs = blobSegmentThresholdGeneral(image,thresholdMethod,...
    filterNoise,filterBackground,minSize,plotRes,saveDir,plotName)
%BLOBSEGMENTTHRESHOLDGENERAL segments blobs in 2D images via various thresholding methods
%
%SYNOPSIS maskBlobs = blobSegmentThreshold(image,minSize,plotRes)
%
%INPUT  image     : 2D image to be segmented.
%       thresholdMethod: 
%                   'otsu' for Otsu.
%                   'rosin' for Rosin.
%                   'minmax' for first minimum after first maximum.
%       filterNoise: 1 to filter out noise, 0 otherwise.
%       filterBackground: 1 to filter out background, 0 otherwise.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       saveDir   : Directory where to save plots.
%                   Only needed if plotRes = 1.
%       plotName  : Image name for plotting.
%                   Only needed if plotRes = 1.
%
%OUTPUT maskBlobs : Mask of blobs. 1 inside blobs, 0 outside.
%
%Khuloud Jaqaman May 2011

%% Output
maskBlobs = [];

%% Input

%check number of input arguments
if nargin < 1
    disp('Please enter at least image to be segmented');
    return
end

%thresholding method
if nargin < 2 || isempty(thresholdMethod)
    thresholdMethod = 'Otsu';
end

%noise filtering
if nargin < 3 || isempty(filterNoise)
    filterNoise = 1;
end

%background filtering
if nargin < 4 || isempty(filterBackground)
    filterBackground = 1;
end

%minimum blob size
if nargin < 5 || isempty(minSize)
    minSize = 20;
end

%plot results
if nargin < 6 || isempty(plotRes)
    plotRes = 0;
end

%mask
mask = ones(size(image));

if ~logical(mask)
    error('Mask must be a logical image.');
end
    
%% Segmentation

%make sure that image is in double format
image = double(image);

%remove noise by filtering image with a Gaussian whose sigma = 1 pixel
if filterNoise
    imageFiltered = filterGauss2D(image,1);
else
    imageFiltered = image;
end

%estimate background by filtering image with a Gaussian whose sigma = 10 pixels
if filterBackground
    imageBackground = filterGauss2D(image,10);
else
    imageBackground = zeros(size(image));
end

%calculate noise-filtered and background-subtracted image
imageFilteredMinusBackground = imageFiltered - imageBackground;

%crop image
imageFilteredMinusBackground = imageFilteredMinusBackground .* mask;

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

%find nonzero values (due to masking)
nzInd = find(imageDilated);

%get minumum and maximum pixel values in image
minSignal = min(imageDilated(nzInd));
maxSignal = max(imageDilated(nzInd));

%normalize nonzero value between 0 and 1
imageDilatedNorm = zeros(size(imageDilated));
imageDilatedNorm(nzInd) = (imageDilated(nzInd) - minSignal) / (maxSignal - minSignal);

%estimate the intensity level to use for thresholding the image
switch thresholdMethod
    case 'otsu'
        try
            level = graythresh(imageDilatedNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (otsu, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'rosin'
        try
            [dummy,level] = cutFirstHistMode(imageDilatedNorm(nzInd),0); %#ok<ASGLU>
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (rosin, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'minmax'
        try
            level = thresholdFluorescenceImage(imageDilatedNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (minmax, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
end

%threshold the image
imageThresholded = im2bw(imageDilatedNorm,level);

%fill holes in thresholded image to make continuous blobs
imageThresholdedFilled = imfill(imageThresholded,'holes');

% go over blobs and remove those with a size smaller that minSize
labels = bwlabel(imageThresholdedFilled);
stats = regionprops(labels, 'Area'); %#ok<MRPBW>
idx = find([stats.Area] > minSize);

%output final blob mask
maskBlobs = ismember(labels, idx);

%% Plotting

if plotRes
    
    %get the blob edges from the final blob mask
    SE = strel('square',3);
    edgesBlobs = imdilate(maskBlobs,SE) - maskBlobs;
    
    %scale the original image to be between 0 and 1
    %actually scale it between the 1st and 99th percentiles
    imageScaled = (image - prctile(image(:),1)) / (prctile(image(:),99) - prctile(image(:),1));
    imageScaled(imageScaled<0) = 0;
    imageScaled(imageScaled>1) = 1;
    
    %plot image
    figHandle1 = figure('Name',[plotName '_segmentation_' ...
        thresholdMethod '_noise' num2str(filterNoise) ...
        '_background' num2str(filterBackground)]);
    hold on
    subplot(1,2,1)
    imshow(imageScaled,[])

    %give the edge pixels a value of zero in the original image
    imageScaled(edgesBlobs==1) = 0;
    
    %construct a 3-layered image to show blob edges on top of
    %original image
    image3Color = repmat(imageScaled,[1 1 3]);
    image3Color(:,:,1) = image3Color(:,:,1) + edgesBlobs;
    
    %plot mask edges
    subplot(1,2,2)
    imshow(image3Color,[]);
%     saveas(figHandle1,fullfile(saveDir,[plotName '_segmentation_' ...
%         thresholdMethod '_noise' num2str(filterNoise) ...
%         '_background' num2str(filterBackground)]),'fig');
    
    %also plot intensity histogram and threshold
    n = histogram(imageDilatedNorm(nzInd),[],0);
    figHandle2 = figure('Name',[plotName '_histogram_' ...
        thresholdMethod '_noise' num2str(filterNoise) ...
        '_background' num2str(filterBackground)]);
    hold on
    histogram(imageDilatedNorm(nzInd),[],0);
    plot(level*[1 1],[0 max(n)],'r','LineWidth',2)
%     saveas(figHandle2,fullfile(saveDir,[plotName '_histogram_' ...
%         thresholdMethod '_noise' num2str(filterNoise) ...
%         '_background' num2str(filterBackground)]),'fig');
%     
%     close all
    
end

%% ~~~ the end ~~~

