function [mask,prctileUsed] = refineEdgeWithSteerableFilterSeed(image,mask0,...
    threshParam,gapCloseParam,doPlot,meanBkg,channel2)
%refineEdgeWithSteerableFilterSeed refines cell edge using intensity gradients obtained from a steerable line filter
%
%SYNOPSIS [mask,prctileUsed] = refineEdgeWithSteerableFilterSeed(image,mask0,...
%    threshParam,gapCloseParam,doPlot,meanBkg,channel2)
%
%INPUT  image        : Image to be segmented.
%       mask0        : Original mask to be refined.
%       threshParam  : Structure with parameters for edge candidate thresholding:
%           .filterSigma    : Standard deviation for steerable filtering.
%                             Optional. Default: 1.5.
%           .prctile        : Percentile for thresholding.
%                             Optional. Default: [95 90 85 80].
%           .bandHalfWidth  : Half-width of band around original edge to
%                             look into for edge refinement.
%                             Input -1 to use whole image instead of a band.
%                             Optional. Default: 50 pixels.
%       gapCloseParam: Structure with parameters for edge gap closing:
%           .maxEdgePairDist: Maximum distance between edge segment pairs.
%                             Optional. Default: 5 pixels.
%           .factorContr    : Contribution of each factor to the edge gap
%                             closing cost. 4 entries for the factors:
%                             (1) distance,
%                             (2) angle between gradients,
%                             (3) angle between gradient and perpendicular
%                                 to centroid-centroid distance,
%                             (4) "edginess" score,
%                             Optional. Default: ones(1,4).
%       doPlot       : 1 to plot masks in the end, 2 to also show edge progress,
%                      0 to plot nothing. In final plot, refined masks
%                      shown in green, original masks shown in blue.
%                      Optional. Default: 0.
%       meanBkg      : Mean background intensity close to the cell edge.
%                      Optional. If not input, information not used.
%       channel2     : Image or some derivative thereof from second
%                      channel if available.
%                      Optional. If not input, information not used.
%
%OUTPUT mask         : Mask (1 inside cell, 0 outside).
%       perctileUsed : Percentile used for gradient thresholding. -1
%                      indicates failed segmentation.
%
%Khuloud Jaqaman, November 2011

%% Input

if nargin < 2
    error('refineEdgeWithSteerableFilterSeed: Wrong number of input arguments');
end

%get thresholding parameters, including steerable filter parameters
if nargin < 3 || isempty(threshParam)
    threshParam.filterSigma = 1.5;
    threshParam.prctile = [95 90 85 80];
    threshParam.bandHalfWidth = 50;
else
    if ~isfield(threshParam,'fiterSigma')
        threshParam.filterSigma = 1.5;
    end
    if ~isfield(threshParam,'prctile')
        threshParam.prctile = [95 90 85 80];
    end
    if ~isfield(threshParam,'bandHalfWidth')
        threshParam.bandHalfWidth = 50;
    end
end
filterSigma = threshParam.filterSigma;
threshPrctile = threshParam.prctile;
bandWidth = 2*threshParam.bandHalfWidth + 1;

%get edge gap closing parameters
if nargin < 4 || isempty(gapCloseParam)
    gapCloseParam.maxEdgePairDist = 5;
    gapCloseParam.factorContr = ones(1,4);
else
    if ~isfield(gapCloseParam,'maxEdgePairDist')
        gapCloseParam.maxEdgePairDist = 5;
    end
    if ~isfield(gapCloseParam,'factorContr')
        gapCloseParam.factorContr = ones(1,4);
    end
end
maxEdgePairDist = gapCloseParam.maxEdgePairDist;
contr = gapCloseParam.factorContr;

%check whether/what to plot
if nargin < 5 || isempty(doPlot)
    doPlot = 0;
end

%get background information
if nargin < 6 || isempty(meanBkg)
    meanBkg = [];
end

%get channel 2 information
if nargin < 7 || isempty(channel2)
    channel2 = [];
end

%% Pre-processing

%convert to double
image = double(image);
mask0 = double(mask0);

%get image size
imSize = size(image);

%get minimum mask size for segmentation assessment
minMaskSize = 0.5*sum(mask0(:));

%get the edges at image boundary of original mask
%and determine whether mask should touch image boundary
mask0Bound = mask0;
mask0Bound(2:end-1,2:end-1) = 0;
if sum(mask0(1,:)) > 0
    mask0Bound(1,3:end-2) = 1;
end
if sum(mask0(end,:)) > 0
    mask0Bound(end,3:end-2) = 1;
end
if sum(mask0(:,1)) > 0
    mask0Bound(3:end-2,1) = 1;
end
if sum(mask0(:,end)) > 0
    mask0Bound(3:end-2,end) = 1;
end

%% Edge candidates and scores

%run steerable filter to enhance edges
[res,theta] = steerableDetector(image,3,filterSigma);

%get the non-maximum suppression image to get edge candidates
nms0 = nonMaximumSuppression(res,theta);

%break junctions in nms image to get clean edge segments
nmsMaskThin = double(bwmorph(nms0~=0, 'thin'));
nn2 = (imfilter(nmsMaskThin,ones(3),'same')-1) .* nmsMaskThin;
junctionMatrix = nn2>2;
nms = nms0 .*nmsMaskThin .* ~junctionMatrix;
% nmsJunctions = nms0 .* nmsMaskThin .* junctionMatrix;

%keep only a certain band of pixels around initial mask
if bandWidth > 0
    goodArea = imdilate(mask0,strel('square',3'))-mask0;
    goodArea = imdilate(goodArea,strel('square',bandWidth));
    nms = nms .* goodArea;
end

%extract some connected component properties in the nms image
nmsBW = (nms > 0);
CC = bwconncomp(nmsBW,8);
numL = CC.NumObjects;
candPixLinAll = CC.PixelIdxList';
[candLengthAll,candMeanGradientAll,candMeanThetaAll,candStdThetaAll,...
    candIntensityAll,candIntSlabUpAll,candIntSlabDnAll,channel2GradAll]...
    = deal(NaN(numL,1));
[candPixAll,candPixSlabUpAll,candPixSlabDnAll] = deal(cell(numL,1));

if numL == 0
    mask = zeros(imSize);
    prctileUsed = -1;
    return
end

for iL = 1 : numL
    
    %pixels of edge itself
    pixL = candPixLinAll{iL};
    [pixLy,pixLx] = ind2sub(imSize,candPixLinAll{iL});
    candPixAll{iL} = [pixLx pixLy]; %pixels
    
    %gradient information
    candMeanGradientAll(iL) = mean(nms(pixL)); %mean gradient at edge
    thetaVal1 = theta(pixL);
    thetaVal2 = thetaVal1;
    thetaVal2(thetaVal2<0) = 2*pi + thetaVal2(thetaVal2<0);
    meanTheta1 = mean(thetaVal1);
    meanTheta2 = mean(thetaVal2);
    varTheta1 = var(thetaVal1);
    varTheta2 = var(thetaVal2);
    if varTheta1 <= varTheta2
        candMeanThetaAll(iL) = meanTheta1; %mean angle at edge
        candStdThetaAll(iL) = sqrt(varTheta1); %std of angle at edge
    else
        if meanTheta2 > pi
            meanTheta2 = meanTheta2 - 2*pi;
        end
        candMeanThetaAll(iL) = meanTheta2; %mean angle
        candStdThetaAll(iL) = sqrt(varTheta2); %std of angle
    end
    
    %pixels of band around edge, up and down the gradient
    if candMeanThetaAll(iL) >= -pi/8 && candMeanThetaAll(iL) <= pi/8 %approximate as 0
        thetaTmp =  0;
    elseif candMeanThetaAll(iL) > pi/8 && candMeanThetaAll(iL) <= 3*pi/8  %approximate as pi/4
        thetaTmp =  pi/4;
    elseif candMeanThetaAll(iL) < -pi/8 && candMeanThetaAll(iL) >= -3*pi/8 %approximate as -pi/4
        thetaTmp = -pi/4;
    elseif candMeanThetaAll(iL) > 3*pi/8 && candMeanThetaAll(iL) <= 5*pi/8 %approximate as pi/2
        thetaTmp =  pi/2;
    elseif candMeanThetaAll(iL) < 3*pi/8 && candMeanThetaAll(iL) >= -5*pi/8 %approximate as -pi/2
        thetaTmp = -pi/2;
    elseif candMeanThetaAll(iL) > 5*pi/8 && candMeanThetaAll(iL) <= 7*pi/8 %approximate as 3pi/4
        thetaTmp =  3*pi/4;
    elseif candMeanThetaAll(iL) < -5*pi/8 && candMeanThetaAll(iL) >= -7*pi/8 %approximate as -3pi/4
        thetaTmp = -3*pi/4;
    else %approximate as pi=-pi
        thetaTmp = pi;
    end
    upX = round(cos(thetaTmp));
    upY = round(sin(thetaTmp));
    
    pixSlabUpX = [pixLx+1*upX; pixLx+2*upX; pixLx+3*upX; pixLx+4*upX; pixLx+5*upX; pixLx+6*upX];
    pixSlabUpY = [pixLy+1*upY; pixLy+2*upY; pixLy+3*upY; pixLy+4*upY; pixLy+5*upY; pixLy+6*upY];
    indx = find(pixSlabUpX<1|pixSlabUpX>imSize(2)|pixSlabUpY<1|pixSlabUpY>imSize(1));
    pixSlabUpX(indx) = [];
    pixSlabUpY(indx) = [];
    candPixSlabUpAll{iL} = [pixSlabUpX pixSlabUpY];
    pixSlabUp = sub2ind(imSize,pixSlabUpY,pixSlabUpX);
    
    pixSlabDnX = [pixLx-1*upX; pixLx-2*upX; pixLx-3*upX; pixLx-4*upX; pixLx-5*upX; pixLx-6*upX];
    pixSlabDnY = [pixLy-1*upY; pixLy-2*upY; pixLy-3*upY; pixLy-4*upY; pixLy-5*upY; pixLy-6*upY];
    indx = find(pixSlabDnX<1|pixSlabDnX>imSize(2)|pixSlabDnY<1|pixSlabDnY>imSize(1));
    pixSlabDnX(indx) = [];
    pixSlabDnY(indx) = [];
    candPixSlabDnAll{iL} = [pixSlabDnX pixSlabDnY];
    pixSlabDn = sub2ind(imSize,pixSlabDnY,pixSlabDnX);
    
    %edge length
    candLengthAll(iL) = length(pixL); %length
    
    %intensity of edge and bands around it
    candIntensityAll(iL) = mean(image(pixL)); %mean intensity at edge
    candIntSlabUpAll(iL) = mean(image(pixSlabUp)); %mean intensity up the gradient
    candIntSlabDnAll(iL) = mean(image(pixSlabDn)); %mean intensity down the gradient
    
    %channel 2 gradient in bands around edge
    if ~isempty(channel2)
        channel2Up = sum(channel2(pixSlabUp));
        channel2Dn = sum(channel2(pixSlabDn));
        channel2GradAll(iL) = channel2Up-channel2Dn;
    end
    
end

%calculate a score for each segment, reflecting how "edgy" it is
if isempty(meanBkg)
    gradOverIntAll = candMeanGradientAll./candIntSlabDnAll;
else
    gradOverIntAll = candMeanGradientAll./max(meanBkg,candIntSlabDnAll);
end
candScoreAll = (candLengthAll.^0.2)/(max(candLengthAll).^0.2) ...
    + gradOverIntAll/max(gradOverIntAll) ...
    - candStdThetaAll/max(candStdThetaAll);
if ~isempty(meanBkg)
    candScoreAll = candScoreAll - abs(candIntSlabDnAll-meanBkg)/max(abs(candIntSlabDnAll-meanBkg));
end
if ~isempty(channel2)
    candScoreAll = candScoreAll + channel2GradAll/max(channel2GradAll);
end

for iPrctile = 1 : length(threshPrctile)
    
    %% Thresholding
    
    %determine threshold for edge segmentation
    cutThreshEdge = prctile(candScoreAll,threshPrctile(iPrctile));
    
    %keep only edge segments above the threshold
    indxEdgeKeep = find(candScoreAll>=cutThreshEdge);
    numEdgeKeep = length(indxEdgeKeep);
    
    %keep information of surviving edges
    candPix = candPixAll(indxEdgeKeep);
    candPixLin = candPixLinAll(indxEdgeKeep);
    candMeanTheta = candMeanThetaAll(indxEdgeKeep);
    candScore = candScoreAll(indxEdgeKeep);
    candIntSlabDn = candIntSlabDnAll(indxEdgeKeep);
    channel2Grad = channel2GradAll(indxEdgeKeep);
    
    %get maximum edge score for later use
    maxScore = max(candScore);
    
    %% Edge linking and gap closing
    
    %Gap closing cost - START
    
    %distance between all surviving edge segment pairs
    %distances = nearest pixel-pixel distance
    %also angle of the nearest-nearest vector
    [edgePairDist,edgePairAnglePerp] = deal(NaN(numEdgeKeep));
    edgePairDistPix = zeros(numEdgeKeep,numEdgeKeep,2,2);
    for iEdge = 1 : numEdgeKeep
        
        %get pixels belonging to iEdge
        pixLIxy = candPix{iEdge};
        
        for jEdge = iEdge + 1 : numEdgeKeep
            
            %get pixels belonging to jEdge
            pixLJxy = candPix{jEdge};
            
            %distance between the two closest pixels
            %             [idx,dist] = KDTreeClosestPoint(pixLJxy,pixLIxy);
            %             kd = KDTree(pixLJxy);
            %             [idx,dist] = nn(kd,pixLIxy);
            [idx,dist] = knnsearch(pixLJxy,pixLIxy);
            
            [dist,indxSort] = sort(dist);
            idx2 = indxSort(1);
            edgePairDist(iEdge,jEdge) = dist(1);            
 
            edgePairDistPix(iEdge,jEdge,:,1) = pixLIxy(idx2(1),:); %image coord
            edgePairDistPix(iEdge,jEdge,:,2) = pixLJxy(idx(idx2(1)),:); %image coord
            
            %angle of vector connecting the two closest pixels
            tmp = pixLJxy(idx(idx2(1)),:) - pixLIxy(idx2(1),:);
            edgePairAnglePerp(iEdge,jEdge) = abs(atan(tmp(2)/tmp(1)));
            
            %symmetrize
            edgePairDistPix(jEdge,iEdge,:,:) = edgePairDistPix(iEdge,jEdge,:,:);
            edgePairDist(jEdge,iEdge) = edgePairDist(iEdge,jEdge);
            edgePairAnglePerp(jEdge,iEdge) = edgePairAnglePerp(iEdge,jEdge);
            
        end
    end
    %angle of perpendicular to vector connecting the two closest pixels
    edgePairAnglePerp = pi/2 - edgePairAnglePerp;
    
    %pair-wise difference in mean theta
    dotProd = cos(candMeanTheta)*cos(candMeanTheta)' + sin(candMeanTheta)*sin(candMeanTheta)';
    dotProd(dotProd<=-1) = -1;
    dotProd(dotProd>=1) = 1;
    meanThetaDiff = acos(dotProd);
    
    %angle between connecting vector and mean theta
    dotProd = abs(cos(abs(repmat(candMeanTheta,1,numEdgeKeep))).*cos(edgePairAnglePerp) + sin(abs(repmat(candMeanTheta,1,numEdgeKeep))).*sin(edgePairAnglePerp));
    dotProd(dotProd>=1) = 1;
    edgePairAngleThetaDiff1 = acos(dotProd);
    dotProd = abs(cos(abs(repmat(candMeanTheta',numEdgeKeep,1))).*cos(edgePairAnglePerp) + sin(abs(repmat(candMeanTheta',numEdgeKeep,1))).*sin(edgePairAnglePerp));
    dotProd(dotProd>=1) = 1;
    edgePairAngleThetaDiff2 = acos(dotProd);
    edgePairAngleThetaDiff = cat(3,edgePairAngleThetaDiff1,edgePairAngleThetaDiff2);
    edgePairAngleThetaDiff = max(edgePairAngleThetaDiff,[],3);
    
    %reward for higher "edginess" score
    pairScoreReward = (repmat(candScore,1,numEdgeKeep) + repmat(candScore',numEdgeKeep,1))/2;
    
    %remove improbable links
    edgePairDist(edgePairDist > maxEdgePairDist) = NaN;
    meanThetaDiff(isnan(edgePairDist)) = NaN;
    edgePairAngleThetaDiff(isnan(edgePairDist)) = NaN;
    pairScoreReward(isnan(edgePairDist)) = NaN;
    
    %calculate linking cost for all pairs
    linkCostAllPairs = contr(1)*edgePairDist/maxEdgePairDist ...
        + contr(2)*meanThetaDiff/pi ...
        + contr(3)*edgePairAngleThetaDiff/(pi/2) ...
        - contr(4)*pairScoreReward/maxScore;
    
    %Gap closing cost - END
    
    %find edge segment with highest score
    %this is most likely a true edge segment, so use it as a seed to
    %put segments together and grow the edge
    scoreSorted = sort(candScore,'descend');
    indxLink = find(candScore==scoreSorted(1));
    indxLink = indxLink(1);
    indxSeed = indxLink;
    
    %start constructing an edge image to assess whether segmentation has been
    %achieved
    imageConst = zeros(imSize);
    imageConst(candPixLin{indxLink}) = 1;
    
    %show edge progress
    if doPlot==2
        figure, hold on
        image3 = cat(3,image/max(image(:)),imageConst,zeros(imSize));
        imshow(image3,[]);
        pause(0.1);
    end
    
    %check whether segmentation has been achieved
    [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize);
    
    iCounter = 0;
    
    %loop and add segments until segmentation has been achieved
    while ~segmentOK && ~isempty(indxLink)
        
        iCounter = iCounter + 1;
        
        %get the relevant costs
        linkCost = linkCostAllPairs(:,indxSeed);
        
        %find possible links
        [indxLink,indxSeed2Link] = find(~isnan(linkCost));
        linkCost = linkCost(sub2ind(size(linkCost),indxLink,indxSeed2Link));
        
        if ~isempty(linkCost)
            
            %take the link with the minimum cost
            indxMin = find(linkCost==min(linkCost));
            indxMin = indxMin(1);
            indxLink = indxLink(indxMin);
            indxSeed2Link = indxSeed(indxSeed2Link(indxMin));
            
            %update list of seeds
            indxSeed = [indxSeed; indxLink]; %#ok<AGROW>
            
            %remove segment pair from possible links
            linkCostAllPairs(indxLink,indxSeed2Link) = NaN;
            
            %add edge segment to edge image
            imageConst(candPixLin{indxLink}) = 1;
            
            %make a link between this edge segment and the existing edge
            %segment that it is linked to
            pixLinkX = squeeze(edgePairDistPix(indxLink,indxSeed2Link,1,:));
            pixLinkY = squeeze(edgePairDistPix(indxLink,indxSeed2Link,2,:));
            pixLinkDist = sqrt(diff(pixLinkX)^2+diff(pixLinkY)^2);
            numSteps = round(pixLinkDist*5);
            pixLinkStepsY = round(linspace(pixLinkX(1),pixLinkX(2),numSteps)); %convert to matrix coord
            pixLinkStepsX = round(linspace(pixLinkY(1),pixLinkY(2),numSteps));
            pixLinkInd = sub2ind(imSize,pixLinkStepsX,pixLinkStepsY);
            imageConst(pixLinkInd) = 1;
            
            %show edge progress
            if doPlot == 2
                image3 = cat(3,image/max(image(:)),imageConst,zeros(imSize));
                imshow(image3,[]);
                pause(0.1);
            end
            
            %check whether full segmentation has been achieved
            if mod(iCounter,10)==0
                [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize);
            end
            
        end
        
    end
    
    %check whether full segmentation has been achieved
    [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize);
    
    %% Final mask
    
    if segmentOK
                
        removeSeed = zeros(numEdgeKeep,1);
        for iSeed = indxSeed'
            
            %copy constructed image to a temporary image
            imageTmp = imageConst;
            
            %remove current edge from temporary image
            imageTmp(candPixLin{iSeed}) = 0;
            
            %check whether segmentation is still OK
            segmentStillOK = checkEdgeComplete(imageTmp,mask0Bound,minMaskSize);
            
            %if OK, indicate that this seed could be removed
            if segmentStillOK
                removeSeed(iSeed) = 1;
            end
            
        end
        indxRemove = find(removeSeed);
        
        %design a metric to sort seeds to be removed
        %this metric goes up as intensity down the gradient goes down
        %it also goes up as number of particles contained in mask goes up
        seedMetric = zeros(size(indxRemove));
        if ~isempty(meanBkg)
            seedMetric = - abs(candIntSlabDn(indxRemove)-meanBkg)/max(abs(candIntSlabDn(indxRemove)-meanBkg));
        end
        if ~isempty(channel2)
            seedMetric = seedMetric + channel2Grad(indxRemove)/max(channel2Grad(indxRemove));
        end
        if all(seedMetric==0)
            seedMetric = candScore(indxRemove);
        end
        
        %order edges in increasing order of this metric
        [~,indxSort] = sort(seedMetric);
        indxRemove = indxRemove(indxSort);
        
        %now go over these seeds and remove seeds that do not break the
        %mask
        for iSeed = indxRemove'
            
            %copy constructed image to a temporary image
            imageTmp = imageConst;
            
            %remove current edge from temporary image
            imageTmp(candPixLin{iSeed}) = 0;
            
            %check whether segmentation is still OK
            [segmentStillOK,imageTmp] = checkEdgeComplete(imageTmp,mask0Bound,minMaskSize);
            
            if segmentStillOK
                
                %if so, delete edge from the constructed image
                imageConst = imageTmp;
                
                %show edge progress
                if doPlot == 2
                    image3 = cat(3,image/max(image(:)),imageConst,zeros(imSize));
                    imshow(image3,[]);
                    pause(0.1);
                end
                
            end
            
        end

        %now make final mask
        edgeMask = imageConst | mask0Bound;
        mask = imfill(edgeMask,'holes');
        
        %get rid of spikes in mask by doing an opening
        SE = strel('disk',3,0);
        mask = imopen(mask,SE);
        
        %also get rid of holes by doing a closure
        SE = strel('disk',3,0);
        mask = imclose(mask,SE);
        
        %if there is more than one connected componenent, retain the
        %largest one
        stats = regionprops(mask,'Area');
        maskArea = vertcat(stats.Area);
        numMasks = length(maskArea);
        if numMasks > 1
            indxKeep = find(maskArea==max(maskArea));
            indxRemove = setdiff(1:numMasks,indxKeep);
            L = bwlabel(mask);
            for iRem = indxRemove
                mask(L==iRem) = 0;
            end
        end
        
        %store percentile used for gradient thresholding
        prctileUsed = threshPrctile(iPrctile);
        
        %exit for loop
        break
        
    end
    
end %(for iPrctile = 1 : length(threshPrctile))

%return empty mask if segmentation is not achieved
if ~segmentOK
    mask = mask0;
    prctileUsed = -1;
end
        
%% Display

%show original mask and modified mask
if doPlot >= 1
    SE = strel('square',3);
    maskEdge0 = mask0 - imerode(mask0,SE);
    maskEdge = mask - imerode(mask,SE);
    imageNorm = (image-min(image(:))) / (max(image(:))-min(image(:)));
    image3 = imageNorm;
    image3(:,:,2) = maskEdge;
    image3(:,:,3) = maskEdge0;
    imtool(image3,[]);
end

%% ~~~ the end ~~~

function [segmentOK,imageConst,mask] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize)

%check whether segmentation has been achieved
edgeMask = imageConst | mask0Bound;
mask = imfill(edgeMask,'holes');
maxArea = sum(mask(:));
if maxArea > minMaskSize
    segmentOK = 1;
else
    segmentOK = 0;
end


