function [mask,prctileUsed] = refineEdgeWithSteerableFilter(mask0,image,...
    threshParam,gapCloseParam,doPlot)
%REFINEEDGEWITHSTEERABLEFILTER refines cell edge using intensity gradients obtained from a steerable line filter
%
%SYNOPSIS [mask,segmentOK] = refineEdgeWithSteerableFilter(mask0,image,...
%    threshParam,gapCloseParam,doPlot)
%
%INPUT  mask0        : Original mask to be refined.
%       image        : Image to be segmented.
%       threshParam  : Structure with parameters for gradient thresholding:
%           .filterSigma    : Standard deviation for filtering.
%                             Optional. Default: 1.5.
%           .gradPrctile    : Gradient percentile for thresholding.
%                             Optional. Default: [95 90 85 80].
%       gapCloseParam: Structure with parameters for edge gap closing:
%           .maxEdgePairDist: Maximum distance between edge segment pair.
%                             Optional. Default: 5 pixels.
%           .maxThetaDiff   : Maximum angle between gradients of edge
%                             segment pair.
%                             Optional. Default: pi (i.e. full range).
%           .maxC2cAngleThetaDiff: Maximum angle between edge gradient and
%                             perpendicular to centroid-centroid vector.
%                             Optional. Default: pi/2 (i.e. full range).
%           .factorContr    : Contribution of each factor to the edge gap
%                             closing cost. 6 entries for the factors:
%                             (1) distance,
%                             (2) angle between gradients,
%                             (3) angle between gradient and perpendicular
%                                 to centroid-centroid distance,
%                             (4) "edginess" score,
%                             (5) intensity,
%                             (6) lack of asymmetry cost.
%                             Optional. Default: [1 1 1 1 1 1].
%           .edgeType       : Flag indicating edge type:
%                             0 = open edge, i.e. image is of part of a
%                             cell and edge touches image boundary.
%                             1 = closed edge, i.e. image is of whole cell
%                             and segmentation requires finding a closed
%                             contour.
%                             2 = closed edge(s), but of potentially more
%                             than one cell. TO BE IMPLEMENTED IN THE
%                             FUTURE IF NEEDED.
%                             Optional. Default: 0.
%           .fracImageCell  : Fraction of image covered by cell. This
%                             number does not have to be accurate, just
%                             some minimum value to help asses whether
%                             segmentation has been achieved.
%                             Optional. Default: 0.25.
%       doPlot       : 1 to plot masks in the end, 2 to also show edge progress,
%                      0 to plot nothing. In final plot, refined masks
%                      shown in green, original masks shown in blue.
%                      Optional. Default: 0.
%
%OUTPUT mask         : Mask (1 inside cell, 0 outside).
%       perctileUsed : Percentile used for gradient thresholding. -1
%                      indicates failed segmentation.
%
%Khuloud Jaqaman, November 2011

%% Input

if nargin < 2
    error('refineEdgeWithSteerableFilter: Wrong number of input arguments');
end

%get thresholding parameters, including steerable filter parameters
if nargin < 3 || isempty(threshParam)
    threshParam.filterSigma = 1.5;
    threshParam.gradPrctile = [95 90 85 80];
else
    if ~isfield(threshParam,'fiterSigma')
        threshParam.filterSigma = 1.5;
    end
    if ~isfield(threshParam,'gradPrctile')
        threshParam.gradPrctile = [95 90 85 80];
    end
end
filterSigma = threshParam.filterSigma;
gradPrctile = threshParam.gradPrctile;

%get edge gap closing parameters
if nargin < 4 || isempty(gapCloseParam)
    gapCloseParam.maxEdgePairDist = 5;
    gapCloseParam.maxThetaDiff = pi;
    gapCloseParam.maxC2cAngleThetaDiff = pi/2;
    gapCloseParam.factorContr = ones(1,6);
    gapCloseParam.edgeType = 0;
    gapCloseParam.fracImageCell = 0.25;
else
    if ~isfield(gapCloseParam,'maxEdgePairDist')
        gapCloseParam.maxEdgePairDist = 5;
    end
    if ~isfield(gapCloseParam,'maxThetaDiff')
        gapCloseParam.maxThetaDiff = pi;
    end
    if ~isfield(gapCloseParam,'maxC2cAngleThetaDiff')
        gapCloseParam.maxC2cAngleThetaDiff = pi/2;
    end
    if ~isfield(gapCloseParam,'factorContr')
        gapCloseParam.factorContr = ones(1,6);
    end
    if ~isfield(gapCloseParam,'edgeType')
        gapCloseParam.edgeType = 0;
    end
    if ~isfield(gapCloseParam,'fracImageCell')
        gapCloseParam.fracImageCell = 0.25;
    end
end
maxEdgePairDist = gapCloseParam.maxEdgePairDist;
maxThetaDiff = gapCloseParam.maxThetaDiff;
maxC2cAngleThetaDiff = gapCloseParam.maxC2cAngleThetaDiff;
contr = gapCloseParam.factorContr;
edgeType = gapCloseParam.edgeType;
fracImageCell = gapCloseParam.fracImageCell;

%make sure that edgeType has one of possible values
if edgeType ~= 0 && edgeType ~= 1
    if edgeType == 2
        error('refineEdgeWithSteerableFilter: Algorithm not developed for edgeType = 2');
    else
        error('refineEdgeWithSteerableFilter: Bad edgeType value');
    end
end

%check whether/what to plot
if nargin < 5 || isempty(doPlot)
    doPlot = 0;
end

%get image size
imSize = size(image);

%get minimum mask size for segmentation assessment
minMaskSize = fracImageCell*prod(imSize);

%% Pre-processing

%run steerable filter to enhance edges
[res,theta] = steerableDetector(image,3,filterSigma);

%run Gaussian filter to smoothen noise
imageF = filterGauss2D(image,filterSigma);

%divide gradient by intensity value to enhance the real cell edge
resOverImage = res./imageF;

%get the non-maximum suppression image to get edge candidates
nmsResOverImage = nonMaximumSuppression(resOverImage,theta);

%get the edges at image boundary of original mask
%needed only in open edge case
mask0Bound = zeros(imSize);
if edgeType == 0
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
    mask0Bound = bwmorph(mask0Bound,'bridge');
end

%extract some connected component properties in the nms image
nmsBW = (nmsResOverImage > 0);
[nmsL,numL] = bwlabel(nmsBW,4);
[candLengthAll,candMeanGradientAll,candMeanThetaAll,candStdThetaAll,...
    candIntensityAll,candMeanGradSlabAll,candIntSlabAll] = deal(NaN(numL,1));
candCentroidAll = NaN(numL,2);
[candPixAll,candPixSlabAll] = deal(cell(numL,1));
SE = strel('square',5);

for iL = 1 : numL
    
    %pixels - edge segment itself
    [pixLy,pixLx] = find(nmsL==iL); %image coord
    candPixAll{iL} = [pixLx pixLy]; %pixels
    pixL = sub2ind(imSize,pixLy,pixLx);
    
    %pixels - slab around edge segment
    imageNew = zeros(imSize);
    imageNew(pixL) = 1;
    imageNew = imdilate(imageNew,SE);
    [pixSlabY,pixSlabX] = find(imageNew); %image coord
    candPixSlabAll{iL} = [pixSlabX pixSlabY]; %pixels
    pixSlab = sub2ind(imSize,pixSlabY,pixSlabX);
    
    %properties - geometrical
    candLengthAll(iL) = length(pixL); %length
    candCentroidAll(iL,:) = [mean(pixLx) mean(pixLy)]; %centroid coord (in image coord)
    
    %properties - gradient
    candMeanGradientAll(iL) = mean(nmsResOverImage(pixL)); %mean gradient at edge
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
    candMeanGradSlabAll(iL) = mean(resOverImage(pixSlab)); %mean gradient in slab
    
    %properties - intensity
    candIntensityAll(iL) = mean(image(pixL)); %mean intensity at edge
    candIntSlabAll(iL) = mean(image(pixSlab)); %mean intensity in slab
    
end

%calculate a score for each segment, reflecting how "edgy" it is
candScoreAll = (candLengthAll.^0.25) .* candMeanGradientAll ./ (1+candStdThetaAll).^2;

for iPrctile = 1 : length(gradPrctile)
    
    %% Thresholding
    
    %determine threshold for edge segmentation
    %     cutThreshEdge = prctile(candMeanGradientAll,gradPrctile(iPrctile));
    cutThreshEdge = prctile(candScoreAll,gradPrctile(iPrctile));
    
    %keep only edge segments above the threshold
    %     indxEdgeKeep = find(candMeanGradientAll>=cutThreshEdge);
    indxEdgeKeep = find(candScoreAll>=cutThreshEdge);
    numEdgeKeep = length(indxEdgeKeep);
    
    %keep information of surviving edges
    candPix = candPixAll(indxEdgeKeep);
    candLength = candLengthAll(indxEdgeKeep);
    candCentroid = candCentroidAll(indxEdgeKeep,:);
    candMeanGradient = candMeanGradientAll(indxEdgeKeep);
    candMeanTheta = candMeanThetaAll(indxEdgeKeep);
    candStdTheta = candStdThetaAll(indxEdgeKeep);
    candIntensity = candIntensityAll(indxEdgeKeep);
    candScore = candScoreAll(indxEdgeKeep);
    
    %get maximum edge intensity and score for later use
    maxInt = max(candIntensityAll);
    maxScore = max(candScoreAll);
    
    %% Edge linking and gap closing
    
    %The different factors contributing to the edge gap closing cost - START
    
    %distance between all surviving edge segment pairs
    %two distances: centroid-centroid & nearest pixel-pixel
    %also angle of the vector connecting the centroids
    [edgePairDist,c2cPairDist,c2cAnglePerp] = deal(NaN(numEdgeKeep));
    edgePairDistPix = zeros(numEdgeKeep,numEdgeKeep,2,2);
    for iEdge = 1 : numEdgeKeep
        
        %get pixels belonging to iEdge
        pixLIxy = candPix{iEdge};
        
        for jEdge = iEdge + 1 : numEdgeKeep
            
            %get pixels belonging to jEdge
            pixLJxy = candPix{jEdge};
            
            %distance between the two closest pixels
            [idx,dist] = KDTreeClosestPoint(pixLJxy,pixLIxy);
            minDistIJ = min(dist);
            idx2 = find(dist==minDistIJ);
            edgePairDist(iEdge,jEdge) = minDistIJ;
            edgePairDistPix(iEdge,jEdge,:,1) = pixLIxy(idx2(1),:); %image coord
            edgePairDistPix(iEdge,jEdge,:,2) = pixLJxy(idx(idx2(1)),:); %image coord
            
            %magnitude and angle of centroid-centroid vector
            tmp = diff(candCentroid([iEdge jEdge],:));
            c2cPairDist(iEdge,jEdge) = norm(tmp);
            c2cAnglePerp(iEdge,jEdge) = abs(atan(tmp(2)/tmp(1)));
            
            %symmetrize
            edgePairDistPix(jEdge,iEdge,:,:) = edgePairDistPix(iEdge,jEdge,:,:);
            edgePairDist(jEdge,iEdge) = edgePairDist(iEdge,jEdge);
            c2cPairDist(jEdge,iEdge) = c2cPairDist(iEdge,jEdge);
            c2cAnglePerp(jEdge,iEdge) = c2cAnglePerp(iEdge,jEdge);
            
        end
    end
    %angle of perpendicular to centroid-centroid vector
    c2cAnglePerp = pi/2 - c2cAnglePerp;
    
    %pair-wise differences in mean gradient and mean theta
    meanGradDiff = abs(repmat(candMeanGradient,1,numEdgeKeep)-repmat(candMeanGradient',numEdgeKeep,1));
    meanGradDiff(isnan(edgePairDist)) = NaN; %#ok<NASGU>
    dotProd = cos(candMeanTheta)*cos(candMeanTheta)' + sin(candMeanTheta)*sin(candMeanTheta)';
    dotProd(dotProd<=-1) = -1;
    dotProd(dotProd>=1) = 1;
    meanThetaDiff = acos(dotProd);
    meanThetaDiff(isnan(edgePairDist)) = NaN;
    
    %angle between centroid-centroid vector and mean theta
    dotProd = abs(cos(abs(repmat(candMeanTheta,1,numEdgeKeep))).*cos(c2cAnglePerp) + sin(abs(repmat(candMeanTheta,1,numEdgeKeep))).*sin(c2cAnglePerp));
    dotProd(dotProd>=1) = 1;
    c2cAngleThetaDiff1 = acos(dotProd);
    dotProd = abs(cos(abs(repmat(candMeanTheta',numEdgeKeep,1))).*cos(c2cAnglePerp) + sin(abs(repmat(candMeanTheta',numEdgeKeep,1))).*sin(c2cAnglePerp));
    dotProd(dotProd>=1) = 1;
    c2cAngleThetaDiff2 = acos(dotProd);
    c2cAngleThetaDiff = cat(3,c2cAngleThetaDiff1,c2cAngleThetaDiff2);
    c2cAngleThetaDiff = max(c2cAngleThetaDiff,[],3);
    
    %reward for higher "edginess" score
    pairScoreReward = (repmat(candScore,1,numEdgeKeep) + repmat(candScore',numEdgeKeep,1))/2;
    
    %reward for higher intensity
    pairIntReward = (repmat(candIntensity,1,numEdgeKeep) + repmat(candIntensity',numEdgeKeep,1))/2;
    
    %links counter to penalize for too many links
    %start at 1 to prevent dividing by zero
    linksCounter = ones(numEdgeKeep);
    
    %The different factors contributing to the edge gap closing cost - END
    
    %remove improbable links
    edgePairDist(edgePairDist > maxEdgePairDist) = NaN;
    meanThetaDiff(meanThetaDiff > maxThetaDiff) = NaN;
    c2cAngleThetaDiff(c2cAngleThetaDiff > maxC2cAngleThetaDiff) = NaN;
    
    %find edge segment with highest score
    %this is most likely a true edge segment, so use it as a seed to
    %put segments together and grow the edge
    if edgeType < 2
        scoreSorted = sort(candScore,'descend');
        indxLink = find(candScore==scoreSorted(1));
        indxLink = indxLink(1);
        indxSeed = indxLink;
    else
        %ADD SOMETHING HERE FOR edgeType = 2
    end
    
    %remove seed from list of segments
    edgePairDist(indxLink,:) = NaN;
    c2cPairDist(indxLink,:) = NaN;
    c2cAngleThetaDiff(indxLink,:) = NaN;
    meanThetaDiff(indxLink,:) = NaN;
    pairScoreReward(indxLink,:) = NaN;
    pairIntReward(indxLink,:) = NaN;
    
    %start constructing an edge image to assess whether segmentation has been
    %achieved
    imageConst = zeros(imSize);
    imageConst(nmsL==indxEdgeKeep(indxLink)) = 1;
    
    %show edge progress
    if doPlot==2
        figure, hold on
        image3 = cat(3,image/max(image(:)),imageConst,zeros(imSize));
        imshow(image3,[]);
        pause(0.1);
    end
    
    %check whether segmentation has been achieved
    %THIS FUNCTION SHOULD BE MODIFIED IN CASE OF WHOLE OR MULTIPLE CELLS
    [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize,edgeType);
    
    %loop and add segments until segmentation has been achieved
    while ~segmentOK && ~isempty(indxLink)
        
        %get the different factors contributing to the costs of potential links
        seedDist         = edgePairDist(:,indxSeed);
        seedC2CPairDist  = c2cPairDist(:,indxSeed); %#ok<NASGU>
        seedC2CAngleDiff = c2cAngleThetaDiff(:,indxSeed);
        seedThetaDiff    = meanThetaDiff(:,indxSeed);
        seedScoreReward  = pairScoreReward(:,indxSeed);
        seedIntReward    = pairIntReward(:,indxSeed);
        includeOrNaN = seedDist;
        includeOrNaN(~isnan(includeOrNaN)) = 1;
        
        %calculate one additional factor, reflecting how asymmetric the
        %distribution of pixels will be when each of the edge segments is
        %added
        %the idea is that the edge contour should be highly asymmetric
        asymmetryCost = NaN(size(seedDist));
        indxInclude = find(max(includeOrNaN,[],2)==1);
        for iPotLink = indxInclude'
            
            %get the edge pixels if this segment is to be added
            newEdgePix = vertcat(candPix{[indxSeed; iPotLink]});
            
            %calculate the positional asymmetry
            pixCovMat = cov(newEdgePix);
            eigenVal = eig(pixCovMat);
            asymmetryCost(iPotLink,:) = min(eigenVal)/max(eigenVal);
            
        end
        
        %calculate the cost
        linkCost = contr(1)*seedDist/maxEdgePairDist ...
            + contr(2)*seedThetaDiff/maxThetaDiff ...
            + contr(3)*seedC2CAngleDiff/maxC2cAngleThetaDiff ...
            - contr(4)*seedScoreReward/maxScore ...
            - contr(5)*seedIntReward/maxInt ...
            + contr(6)*asymmetryCost/max(asymmetryCost(:));
        
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
            
            %remove added segment from possible links
            edgePairDist(indxLink,:) = NaN;
            c2cPairDist(indxLink,:) = NaN;
            c2cAngleThetaDiff(indxLink,:) = NaN;
            meanThetaDiff(indxLink,:) = NaN;
            pairScoreReward(indxLink,:) = NaN;
            pairIntReward(indxLink,:) = NaN;
            
            %add 1 to counter of seed that got linked to
            linksCounter(:,indxSeed2Link) = linksCounter(:,indxSeed2Link) + 1;
            
            %add edge segment to edge image
            imageConst(nmsL==indxEdgeKeep(indxLink)) = 1;
            
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
            %THIS FUNCTION SHOULD BE MODIFIED IN CASE OF WHOLE OR MULTIPLE CELLS
            [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize,edgeType);
            
        end
        
    end
    
    %% Final mask
    
    if segmentOK
        
        %order edges in increasing order of their score
        seedScore = candScore(indxSeed);
        [~,indxSort] = sort(seedScore);
        indxSeed = indxSeed(indxSort);
        
        %go over all edges in order of increasing score
        for iSeed = indxSeed'
            
            %copy constructed image to a temporary image
            imageTmp = imageConst;
            
            %remove current edge from temporary image
            imageTmp(nmsL==indxEdgeKeep(iSeed)) = 0;
            
            %check whether segmentation is still OK
            %THIS FUNCTION SHOULD BE MODIFIED IN CASE OF WHOLE OR MULTIPLE CELLS
            [segmentStillOK,imageTmp] = checkEdgeComplete(imageTmp,mask0Bound,minMaskSize,edgeType);
            
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
        SE = strel('disk',5,0);
        mask = imopen(mask,SE);
        
        %also get rid os whole by doing a closure
        SE = strel('disk',5,0);
        mask = imclose(mask,SE);
        
        %store percentile used for gradient thresholding
        prctileUsed = gradPrctile(iPrctile);
        
        %exit for loop
        break
        
    end
    
end %(for iPrctile = 1 : length(gradPrctile))

%return empty mask if segmentation is not achieved
if ~segmentOK
    mask = zeros(imSize);
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

function [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize,edgeType)

switch edgeType
    
    case 0
        
        %check whether segmentation has been achieved
        edgeMask = imageConst | mask0Bound;
        mask = imfill(edgeMask,'holes');
        maxArea = sum(mask(:));
        if maxArea > minMaskSize
            segmentOK = 1;
        else
            segmentOK = 0;
        end
        
    case 1
        
        %POTENTIALLY ADD SOMETHING HERE TO BRIDGE THAT FINAL GAP IN
        %imageConst (THUS RETURN imageConst AFTER MODIFICATION)
        
        %check whether segmentation has been achieved
        edgeMask = imageConst | mask0Bound;
        mask = imfill(edgeMask,'holes');
        maxArea = sum(mask(:));
        if maxArea > minMaskSize
            segmentOK = 1;
        else
            segmentOK = 0;
        end
        
    case 2
        
        %ADD SOMETHING IN CASE OF MULTIPLE CELLS
        
end


