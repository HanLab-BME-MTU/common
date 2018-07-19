function [labelImg,nLabel] = segmentNucleii(img, dilateNucleii)
%[labelImg] = segmentNucleii(img) sements img into label indentifying each
% nucleus by double thresholding and watershed.
%  input: img - 2D image matrix
% output: labelImg - label matrix (uint16 or uint8) 
%            nLabel - the total number of nucleii
if nargin<2
    dilateNucleii=false;
end
img = double(img);
%% Do background subtraction
%remove noise by filtering image with a Gaussian whose sigma = 1 pixel
imageFiltered = filterGauss2D(img,1);

%estimate background by filtering image with a Gaussian whose sigma = 10 pixels
imageBackground = filterGauss2D(img,10);

%calculate noise-filtered and background-subtracted image
imageFilteredMinusBackground = imageFiltered - imageBackground;
imageDilated = imageFilteredMinusBackground;

%get minumum and maximum pixel values in image
minSignal = min(imageDilated(:));
maxSignal = max(imageDilated(:));
imageDilatedNorm = (imageDilated - minSignal) / (maxSignal - minSignal);

% bwDAPI = im2bw(imgDAPInorm,levelrosin);
% figure, imshow(bwDAPI)
%% Thresholding
% imgnorm = img/max(img(:));
levelOtsu = graythresh(imageDilatedNorm);
[~, level2] = cutFirstHistMode(imageDilatedNorm,0); %Rosin
levelOtsu = 0.99*levelOtsu + 0.01*level2;
levelrosin = 0.33333*levelOtsu + 0.66667*level2;
% bwDAPI = im2bw(imageNormFiltered,levelOtsu);
bwDAPI = im2bw(imageDilatedNorm,levelOtsu);
%% Lable each segments and break the clumps
labelOtsu = bwlabel(bwDAPI);
labelDAPIws = labelOtsu;
areas = regionprops(bwDAPI,'Area');
medianArea = median(arrayfun(@(x) (x.Area),areas));
% Sometimes otsu algorithm can pick up too generous value if there is some
% out-of-focus signal in the image. In that case, we will pass this image
% by returning NaNs for labelImg and nLabel.
maxArea = max(arrayfun(@(x) (x.Area),areas));
if maxArea>1000*medianArea
    disp('Some out-of-focus signal interfered automatic thresholding, leading to huge wrong segmentation. Passing on this cell ...')
    labelImg=NaN;
    nLabel = NaN;
    return
end
% For each segment, check if they can be divided by watershed
nLabel = max(labelOtsu(:));
idxLargeNucleii=find(arrayfun(@(x) x.Area>medianArea,areas));
newLabel = nLabel;
p=0;
nLargeNucleii=length(idxLargeNucleii);
% remove noise by Gaussian filter with sigma of 1
imageDAPIFiltered = filterGauss2D(img,1);

progressText(0,'Splitting clumps: ')
for ii=idxLargeNucleii'
    p=p+1;
    curBW = labelOtsu==ii;
%     curDist = bwdist(~curBW);
%     curDist = -curDist;
%     curDist(~curBW)=-Inf;
%     curDAPI = curBW.*(-imgDAPI);
    curDAPI = curBW.*(-imageDAPIFiltered);
    
    curDAPI(~curBW)=-Inf;
   
%     curWS = watershed(curDist);
    curWS = watershed(curDAPI);
    % 0 is boundary, 1 is background, 2 and above identifies each object
    maxWSLabel = max(curWS(:));
    if maxWSLabel>2
        % Replace labelDAPI with watersheded one
        % empty the corresponding label
        labelDAPIws(curBW)=0;
        for jj=2:maxWSLabel
            curWScomponent = curWS==jj;
            if jj==2
                % Leave existing label as is for WSlabel=2
                labelDAPIws = labelDAPIws+curWScomponent*ii;
            else
                % Make new label above the max label and add it
                newLabel = newLabel+1;
                labelDAPIws = labelDAPIws+curWScomponent*newLabel;
            end
        end
    end
    progressText(p/(nLargeNucleii-1),'Splitting clumps')
end
%% Showing
% figure, imshow(curBW)
% figure, imshow(curDist,[])
% figure, imshow(curWS,[]),colormap jet
% clumpColor = distinguishable_colors(maxWSLabel-1,[0 0 0]);
% clumpColor = [[0 0 0]; clumpColor];
% figure, imshow(curWS,clumpColor)
% watershed for clumps
% figure, imshow(labelOtsu,[]), colormap jet
% figure, imshow(labelDAPIws,[]), colormap jet
%% Work on image dilation for cells bigger than nucleus size
% use labelDAPIws as a seed point to propagate the same IDs to expanded
% version of mask. I can either dilate the existing or use Rosin
% thresholding to get more generous mask. I'll try latter first.
% bwDAPIrosin = im2bw(imageNormFiltered,levelrosin);
bwDAPIrosin = im2bw(imageDilatedNorm,levelrosin);
if dilateNucleii
    dilationRadius = 5; % pixel
    se = strel('disk',dilationRadius);
    bwDAPIrosin = imdilate(bwDAPIrosin,se);
end
% There is some unwanted additional segmentation from lower threshold by
% Rosin. We need to clear those out before processing

% figure, imshow(bwDAPIrosin)
% Now propagate label to bwDAPIrosin
% For each label, calculate distance from the segment in bwDPIrosin. 
% oldNLable = nLabel;
nLabel=newLabel;
distFromAllSeg = zeros(size(labelDAPIws,1),size(labelDAPIws,2),nLabel, 'uint8');
%% Actual assigning
progressText(0,'Assigning ids to dilated mask')
for ii=1:nLabel
    curBW = labelDAPIws==ii;
    curDist = bwdist(curBW);
    distFromAllSeg(:,:,ii) = curDist;
    progressText(ii/(nLabel-1),'Assigning ids to dilated mask')
end
%% Assign the existing pixel an existing label from the labelDAPIws to
% labelDAPIrosin
labelImg = labelDAPIws;
% For signal pixels excluding labelDAPIws, look at the distance and put the
% label that gives the least distance to the pixel.
pixelAdditional = bwDAPIrosin & (~(labelDAPIws>0));
idPixelAdditional = find(pixelAdditional);
tooLongDist = 100;
%% Actual assigning
for k=idPixelAdditional'
    [I,J] = ind2sub(size(bwDAPIrosin),k);
    [minDist,minLabel]=min(distFromAllSeg(I,J,:));
    if minDist<tooLongDist
        labelImg(k)=minLabel;
    end
end
%% Do regionprops and see if there is any double or triple segmentations with the same ID. If then leave the segmentation with largest area
progressText(0,'Removing noise')
for ii=1:nLabel
    curCC = bwconncomp(labelImg==ii);
    if curCC.NumObjects>1
        regionIndNuc = regionprops(curCC,'Area','PixelID');
        [~,indMaxArea]=max(arrayfun(@(x) x.Area,regionIndNuc));
        % Replace minor areas with zeros
        indMinorAreas = setdiff(1:numel(regionIndNuc),indMaxArea);
        for jj=indMinorAreas
            labelImg(regionIndNuc(jj).PixelIdxList)=0;
        end
    end    
    progressText(ii/(nLabel-1),'Removing noise')
end
