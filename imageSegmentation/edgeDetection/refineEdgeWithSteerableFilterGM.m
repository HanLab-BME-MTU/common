function [mask,prctileUsed] = refineEdgeWithSteerableFilterGM(image,mask0,...
    threshParam,gapCloseParam,doPlot,particleInfo,meanBkg)
%refineEdgeWithSteerableFilterGM *OUTDATED* refines cell edge using intensity gradients obtained from a steerable line filter
%
%SYNOPSIS [mask,segmentOK] = refineEdgeWithSteerableFilterGM(mask0,image,...
%    threshParam,gapCloseParam,doPlot)
%
%INPUT  image        : Image to be segmented.
%       mask0        : Original mask to be refined.
%       threshParam  : Structure with parameters for gradient thresholding:
%           .filterSigma    : Standard deviation for filtering.
%                             Optional. Default: 1.5.
%           .gradPrctile    : Gradient percentile for thresholding.
%                             Optional. Default: [95 90 85 80].
%           .bandWidth      : Width of band around original edge to look
%                             into for edge refinement.
%                             Use -1 so as to use whole image instead of a
%                             band.
%                             Optional. Default: 100 pixels (50 on each side).
%       gapCloseParam: Structure with parameters for edge gap closing:
%           .maxEdgePairDist: Maximum distance between edge segment pair.
%                             Optional. Default: 5 pixels.
%           .factorContr    : Contribution of each factor to the edge gap
%                             closing cost. 6 entries for the factors:
%                             (1) distance,
%                             (2) angle between gradients,
%                             (3) angle between gradient and perpendicular
%                                 to centroid-centroid distance,
%                             (4) "edginess" score,
%                             Optional. Default: ones(1,4).
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
%       particleInfo : Detection output with particle positions for
%                      relevant frames.
%                      Optional. If not input, information not used.
%       meanBkg      : Mean background intensity close to the cell edge.
%                      Optional. If not input, information not used.
%
%OUTPUT mask         : Mask (1 inside cell, 0 outside).
%       perctileUsed : Percentile used for gradient thresholding. -1
%                      indicates failed segmentation.
%
%REMARKS NOT UP-TO-DATE -- FOLLOW refineEdgeWithSteerableFilterSeed AND FIX
%        BEFORE USING.
%
%Khuloud Jaqaman, November 2011

%% Input

if nargin < 2
    error('refineEdgeWithSteerableFilterGM: Wrong number of input arguments');
end

%get thresholding parameters, including steerable filter parameters
if nargin < 3 || isempty(threshParam)
    threshParam.filterSigma = 1.5;
    threshParam.gradPrctile = [95 90 85 80];
    threshParam.bandWidth = 100;
else
    if ~isfield(threshParam,'fiterSigma')
        threshParam.filterSigma = 1.5;
    end
    if ~isfield(threshParam,'gradPrctile')
        threshParam.gradPrctile = [95 90 85 80];
    end
    if ~isfield(threshParam,'bandWidth')
        threshParam.bandWidth = 100;
    end
end
filterSigma = threshParam.filterSigma;
gradPrctile = threshParam.gradPrctile;
bandWidth = threshParam.bandWidth;

%get edge gap closing parameters
if nargin < 4 || isempty(gapCloseParam)
    gapCloseParam.maxEdgePairDist = 5;
    gapCloseParam.factorContr = ones(1,4);
    gapCloseParam.edgeType = 0;
    gapCloseParam.fracImageCell = 0.25;
else
    if ~isfield(gapCloseParam,'maxEdgePairDist')
        gapCloseParam.maxEdgePairDist = 5;
    end
    if ~isfield(gapCloseParam,'factorContr')
        gapCloseParam.factorContr = ones(1,4);
    end
    if ~isfield(gapCloseParam,'edgeType')
        gapCloseParam.edgeType = 0;
    end
    if ~isfield(gapCloseParam,'fracImageCell')
        gapCloseParam.fracImageCell = 0.25;
    end
end
maxEdgePairDist = gapCloseParam.maxEdgePairDist;
contr = gapCloseParam.factorContr;
edgeType = gapCloseParam.edgeType;
fracImageCell = gapCloseParam.fracImageCell;

%make sure that edgeType has one of possible values
if edgeType ~= 0 && edgeType ~= 1
    if edgeType == 2
        error('refineEdgeWithSteerableFilterGM: Algorithm not developed for edgeType = 2');
    else
        error('refineEdgeWithSteerableFilterGM: Bad edgeType value');
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

%get the particle positions
if nargin < 6 || isempty(particleInfo)
    linIndxPart = [];
else
    xCoordPart = vertcat(particleInfo.xCoord);
    yCoordPart = vertcat(particleInfo.yCoord);
    linIndxPart = sub2ind(imSize,ceil(yCoordPart(:,1)),ceil(xCoordPart(:,1)));
end

if nargin < 7 || isempty(meanBkg)
    meanBkg = [];
end

%% Pre-processing

%run steerable filter to enhance edges
[res,theta] = steerableDetector(image,3,filterSigma);

% %run Gaussian filter to smoothen noise
% imageF = filterGauss2D(image,filterSigma);

% %divide gradient by intensity value to enhance the real cell edge
% resOverImage = res./imageF;
resOverImage = res;

%get the non-maximum suppression image to get edge candidates
nmsResOverImage0 = nonMaximumSuppression(resOverImage,theta);

%break junctions in nms image to get clean edge segments
nmsMaskThin = double(bwmorph(nmsResOverImage0~=0, 'thin'));
nn = (imfilter(nmsMaskThin,ones(3),'same')-1) .* nmsMaskThin;
junctionMatrix = nn>2;
nmsResOverImage = nmsResOverImage0 .*nmsMaskThin .* ~junctionMatrix;
% nmsJunctions = nmsResOverImage0 .* nmsMaskThin .* junctionMatrix;

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

%keep only a certain band of pixels around initial mask
if bandWidth > 0
    goodArea = imdilate(mask0,strel('square',3'))-mask0;
    goodArea = imdilate(goodArea,strel('square',bandWidth));
    nmsResOverImage = nmsResOverImage .* goodArea;
end

%extract some connected component properties in the nms image
nmsBW = (nmsResOverImage > 0);
CC = bwconncomp(nmsBW,8);
numL = CC.NumObjects;
candPixLinAll = CC.PixelIdxList';
[candLengthAll,candMeanGradientAll,candMeanThetaAll,candStdThetaAll,...
    candIntensityAll,candIntSlabUpAll,candIntSlabDnAll,candParticleGradAll]...
    = deal(NaN(numL,1));
[candPixAll,candPixSlabUpAll,candPixSlabDnAll] = deal(cell(numL,1));
maskUpDn = zeros(imSize);

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
    
    %particle density gradient in bands around edge
    maskUpDn(pixSlabUp) = 1;
    numPartUp = sum(maskUpDn(linIndxPart));
    maskUpDn(pixSlabUp) = 0;
    maskUpDn(pixSlabDn) = 1;
    numPartDn = sum(maskUpDn(linIndxPart));
    maskUpDn(pixSlabDn) = 0;
    tmp = (numPartUp-numPartDn)/mean([numPartUp numPartDn]);
    if isnan(tmp) || isinf(tmp)
        candParticleGradAll(iL) = 0;
    else
        candParticleGradAll(iL) = tmp;
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
if ~isempty(linIndxPart)
    candScoreAll = candScoreAll + candParticleGradAll/max(candParticleGradAll);
end

% candScoreAll = (candLengthAll.^0.25) .* gradOverIntAll ./ (1+candStdThetaAll).^2;
% if ~isempty(linIndxPart)
%     candScoreAll = candScoreAll .* candParticleGradAll;
% end

for iPrctile = 1 : length(gradPrctile)
    
    %% Thresholding
    
    %determine threshold for edge segmentation
    cutThreshEdge = prctile(candScoreAll,gradPrctile(iPrctile));
    
    %keep only edge segments above the threshold
    indxEdgeKeep = find(candScoreAll>=cutThreshEdge);
    numEdgeKeep = length(indxEdgeKeep);
    
    %keep information of surviving edges
    candPix = candPixAll(indxEdgeKeep);
    candPixLin = candPixLinAll(indxEdgeKeep);
    candMeanTheta = candMeanThetaAll(indxEdgeKeep);
    candScore = candScoreAll(indxEdgeKeep);
    candIntSlabDn = candIntSlabDnAll(indxEdgeKeep);
    candParticleGrad = candParticleGradAll(indxEdgeKeep);
    
    %get maximum edge intensity and score for later use
    maxScore = max(candScoreAll);
    
    %% Edge linking and gap closing
    
    %GAP CLOSING COST - START
    
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
            [idx,dist] = KDTreeClosestPoint(pixLJxy,pixLIxy);
            
            %             minDistIJ = min(dist);
            %             idx2 = find(dist==minDistIJ);
            %             edgePairDist(iEdge,jEdge) = minDistIJ;
            
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
    
    for i=1:size(linkCostAllPairs,1)
        linkCostAllPairs(i,i:end) = NaN;
    end
    
    %convert the matrix of linking costs to a list of nodes and edges for
    %weighted graph matching
    [edgePairIdx(:,1),edgePairIdx(:,2)] = find(~isnan(linkCostAllPairs));
    linkReward = -linkCostAllPairs(~isnan(linkCostAllPairs));
    numEdgePairs = length(linkReward);
    
    %rank the links from highest to lowest reward
    [linkReward,indxSort] = sort(linkReward,'descend');
    edgePairIdx = edgePairIdx(indxSort,:);
    
    %define indices of possible links - since this is the very first
    %matching, all links are possible
    indxPossible = (1 : numEdgePairs)';
    
    %GAP CLOSING COST - END
    
    imageConst = zeros(imSize); %this is to check whether segmentation has been achieved
    iCounter = 0;
    indxSeed = [];
    segmentOK = 0;
    
    %loop and add segments until segmentation has been achieved
    while ~segmentOK && (numEdgePairs > 0)
        
        iCounter = iCounter + 1;
        
        %find matches
        edgePairIdxPoss = edgePairIdx(indxPossible,:);
        linkRewardPoss = linkReward(indxPossible);
        edgeMatches = maxWeightedMatching(numEdgeKeep,edgePairIdxPoss,linkRewardPoss);
        indxMatches = indxPossible(edgeMatches);
        
        %         edgeMatches = maxWeightedMatching(numEdgePairs,edgePairIdx,linkReward);
        %         indxMatches = find(edgeMatches);
        
        numMatches = length(indxMatches);
        
        if numMatches > 0
            
            %take only the top 5*(number of cells) matches
            indxMatches = indxMatches(1:min(5*max(edgeType,1),numMatches));
            numMatches = length(indxMatches);
            
            %update list of seeds
            edgePairsMatched = edgePairIdx(indxMatches,:);
            indxSeed = unique([indxSeed; edgePairsMatched(:)]);
            
            %remove matched segments from possible links
            edgePairIdx(indxMatches,:) = [];
            linkReward(indxMatches) = [];
            numEdgePairs = numEdgePairs - numMatches;
            
            %update list of possible links, from now on it has to include
            %segments already in the list of seeds
            indxPossible = [];
            for iSeed = indxSeed'
                indxPossible = [indxPossible; find( edgePairIdx(:,1)==iSeed | edgePairIdx(:,2)==iSeed )]; %#ok<AGROW>
            end
            indxPossible = unique(indxPossible);
            
            %add edge segments to edge image
            imageConst(vertcat(candPixLin{indxSeed})) = 1;
            
            %interpolate between the matched edge segments
            for iMatch = 1 : numMatches
                pixLinkX = squeeze(edgePairDistPix(edgePairsMatched(iMatch,1),edgePairsMatched(iMatch,2),1,:));
                pixLinkY = squeeze(edgePairDistPix(edgePairsMatched(iMatch,1),edgePairsMatched(iMatch,2),2,:));
                pixLinkDist = sqrt(diff(pixLinkX)^2+diff(pixLinkY)^2);
                numSteps = round(pixLinkDist*5);
                pixLinkStepsY = round(linspace(pixLinkX(1),pixLinkX(2),numSteps)); %convert to matrix coord
                pixLinkStepsX = round(linspace(pixLinkY(1),pixLinkY(2),numSteps));
                pixLinkInd = sub2ind(imSize,pixLinkStepsX,pixLinkStepsY);
                imageConst(pixLinkInd) = 1;
            end
            
            %show edge progress
            if doPlot == 2
                image3 = cat(3,image/max(image(:)),imageConst,zeros(imSize));
                imshow(image3,[]);
                pause(0.1);
            end
            
            %check whether full segmentation has been achieved
            %THIS FUNCTION SHOULD BE MODIFIED IN CASE OF WHOLE OR MULTIPLE CELLS
            %         if mod(iCounter,10)==0
            [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize,edgeType);
            %         end
            
        else
            
            continue
            
        end
        
    end
    
    %check whether full segmentation has been achieved
    %THIS FUNCTION SHOULD BE MODIFIED IN CASE OF WHOLE OR MULTIPLE CELLS
    [segmentOK,imageConst] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize,edgeType);
    
    %% Final mask
    
    if segmentOK
        
        removeSeed = zeros(numEdgeKeep,1);
        for iSeed = indxSeed'
            
            %copy constructed image to a temporary image
            imageTmp = imageConst;
            
            %remove current edge from temporary image
            imageTmp(candPixLin{iSeed}) = 0;
            
            %check whether segmentation is still OK
            %THIS FUNCTION SHOULD BE MODIFIED IN CASE OF WHOLE OR MULTIPLE CELLS
            segmentStillOK = checkEdgeComplete(imageTmp,mask0Bound,minMaskSize,edgeType);
            
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
        if ~isempty(linIndxPart)
            seedMetric = seedMetric + candParticleGrad(indxRemove)/max(candParticleGrad(indxRemove));
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
        SE = strel('disk',3,0);
        mask = imopen(mask,SE);
        
        %also get rid os whole by doing a closure
        SE = strel('disk',3,0);
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

function [segmentOK,imageConst,mask] = checkEdgeComplete(imageConst,mask0Bound,minMaskSize,edgeType)

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
