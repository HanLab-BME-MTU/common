function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
     errFlag] = costMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,...
     costMatParam,useLocalDensity,nnDistTracks,probDim,linearMotion,prevCost)
%COSTMATLINEARMOTIONLINK provides a cost matrix for linking features based on competing linear motion models
%
%SYNOPSIS [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
%     errFlag] = costMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,...
%     costMatParam,useLocalDensity,nnDistTracks,probDim,linearMotion,prevCost)
%
%INPUT  movieInfo             : A 2x1 array (corresponding to the 2 frames of 
%                               interest) containing the fields:
%             .allCoord           : x,dx,y,dy,[z,dz] of features collected in one
%                                   matrix.
%             .amp                : Amplitudes of PSFs fitting detected features. 
%                                   1st column for values and 2nd column 
%                                   for standard deviations.
%             .num                : Number of features in each frame.
%             .nnDist             : Distance from each feature to its nearest
%                                   neighbor. Not needed at the moment.
%      kalmanFilterInfoFrame1 : Structure with the following fields:
%             .stateVec           : Kalman filter state vector for each
%                                   feature in 1st frame.
%             .stateCov           : Kalman filter state covariance matrix
%                                   for each feature in 1st frame.
%             .noiseVar           : Variance of state noise for each
%                                   feature in 1st frame.
%      costMatParam           : Structure with fields:
%             .maxSearchRadiusL   : Maximum allowed search radius (in pixels).
%             .minSearchRadiusL   : Minimum allowed search radius (in pixels).
%             .brownStdMultL      : Factor multiplying Brownian
%                                   displacement std to get search radius.
%             .closestDistScaleL  :Scaling factor of nearest neighbor
%                                   distance. Not needed if useLocalDensity = 0;
%             .maxStdMultL        : Maximum value of factor multiplying
%                                   Brownian displacement std to get search
%                                   radius. Not needed if useLocalDensity = 0;
%      useLocalDensity        : Logical variable indicating whether to use
%                               local density in search radius estimation.
%      nnDistTracks           : Nearest neighbor distance of features in
%                               frame 1 given their history.
%      probDim                : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%      linearMotion           : 1 if linear motion is to be considered, 0 
%                               otherwise.
%      prevCost               : Matrix of previous linking costs.
%
%OUTPUT costMat               : Cost matrix.
%       propagationScheme     : Propagation scheme corresponding to each
%                               cost in the cost matrix.
%       kalmanFilterInfoFrame2: Structure with the following fields:
%             .stateVec           : Kalman filter prediction of state
%                                   vector in 2nd frame based on all 3 
%                                   motion models.
%             .stateCov           : Kalman filter prediction of state
%                                   covariance matrix in 2nd frame based on
%                                   all 3 motion models.
%             .obsVec             : Kalman filter prediction of the
%                                   observed variables in 2nd frame based
%                                   on all 3 motion models.
%       nonlinkMarker         : Value indicating that a link is not allowed.
%       errFlag               : 0 if function executes normally, 1 otherwise.
%
%REMARKS Three competing linear motion models: 1, 2 and 3.
%        1: forward drift, 2: backward drift, 3: zero drift (Brownian).
%
%Khuloud Jaqaman, March 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costMat = [];
propagationScheme = [];
kalmanFilterInfoFrame2 = [];
nonlinkMarker = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatLinearMotionLink')
    disp('--costMatLinearMotionLink: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cost calculation parameters
maxSearchRadius = costMatParam.maxSearchRadiusL;
minSearchRadius = costMatParam.minSearchRadiusL;
brownStdMult    = costMatParam.brownStdMultL;
if useLocalDensity
    closestDistScale = costMatParam.closestDistScaleL;
    maxStdMult = costMatParam.maxStdMultL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify number of propagation schemes used
numSchemes = 3;

%calculate vector sizes
vecSize = 2 * probDim;

%construct transition matrices
if linearMotion
    transMat(:,:,1) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(-ones(probDim,1),probDim); %backward drift transition matrix
    transMat(:,:,3) = eye(vecSize); %zero drift transition matrix
else
    transMat(:,:,3) = eye(vecSize) + diag(ones(probDim,1),probDim); %forward drift transition matrix
    transMat(:,:,2) = eye(vecSize) + diag(-ones(probDim,1),probDim); %backward drift transition matrix
    transMat(:,:,1) = eye(vecSize); %zero drift transition matrix
end
%construct observation matrix
observationMat = [eye(probDim) zeros(probDim)]; %observation matrix

%get number of features in the 2 frames
numFeaturesFrame1 = movieInfo(1).num;
numFeaturesFrame2 = movieInfo(2).num;

%reserve memory for "kalmanFilterInfoframe2"
kalmanFilterInfoFrame2 = struct('stateVec',zeros(numFeaturesFrame1,vecSize,numSchemes),...
    'stateCov',zeros(vecSize,vecSize,numFeaturesFrame1,numSchemes),...
    'obsVec',zeros(numFeaturesFrame1,probDim,numSchemes));

%apply Kalman filters to each feature in 1st frame
for iFeature = 1 : numFeaturesFrame1

    %get state vector and its covariance matrix of feature in 1st frame
    stateOld = kalmanFilterInfoFrame1.stateVec(iFeature,:)';
    stateCovOld = kalmanFilterInfoFrame1.stateCov(:,:,iFeature);
    noiseVar = kalmanFilterInfoFrame1.noiseVar(:,:,iFeature);

    %go over all possible propagation schemes
    for iScheme = 1 : numSchemes

        %predict state vector of feature in 2nd frame
        stateVec = transMat(:,:,iScheme)*stateOld;

        %predict state covariance matrix of feature in 2nd frame
        stateCov = transMat(:,:,iScheme)*stateCovOld*transMat(:,:,iScheme)' ...
            + noiseVar;

        %determine observation vector of feature in 2nd frame (i.e. the
        %propagated position of the feature)
        obsVec = observationMat*stateVec;

        %save information in kalmanFilterInfoFrame2
        kalmanFilterInfoFrame2.stateVec(iFeature,:,iScheme) = stateVec';
        kalmanFilterInfoFrame2.stateCov(:,:,iFeature,iScheme) = stateCov;
        kalmanFilterInfoFrame2.obsVec(iFeature,:,iScheme) = obsVec';

    end

end

%get the propagated positions of features in 1st frame based on the three propagation schemes
propagatedPos = kalmanFilterInfoFrame2.obsVec;

%put the coordinates of features in the 2nd frame in one matrix
coord2 = movieInfo(2).allCoord(:,1:2:end);

%calculate the cost matrices for all three propagation schemes
for iScheme = 1 : numSchemes

    %put the propagated x and y coordinates of features from 1st frame in
    %one matrix
    coord1 = propagatedPos(:,:,iScheme);

    %calculate the distances between features
    costMatTmp(:,:,iScheme) = createDistanceMatrix(coord1,coord2);

end

%find the minimum cost for the link between every pair, which also
%determines the best propagation scheme to perform that link
[costMat,propagationScheme] = min(costMatTmp,[],3);

%get the Kalman standard deviation of all features in frame 1
kalmanStd = sqrt(probDim * squeeze(kalmanFilterInfoFrame1.noiseVar(1,1,:)));

%copy brownStdMult into vector
stdMultInd = repmat(brownStdMult,numFeaturesFrame1,1);

%if local density information is used to expand search radius ...
if useLocalDensity

    %divide each feature's nearest neighbor distance/closestDistScale by kalmanStd
    ratioDist2Std = nnDistTracks./kalmanStd/closestDistScale;

    %make ratios larger than maxStdMult equal to maxStdMult
    ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

    %expand search radius multiplication factor if possible
    stdMultInd = max([stdMultInd ratioDist2Std],[],2);
    
end

%get the search radius of each feature in frame 1 and make sure it falls
%within reasonable limits
searchRadius = stdMultInd .* kalmanStd;
searchRadius(searchRadius>maxSearchRadius) = maxSearchRadius;
searchRadius(searchRadius<minSearchRadius) = minSearchRadius;

%replicate the search radius to compare to cost matrix
searchRadius = repmat(searchRadius,1,numFeaturesFrame2);

%assign NaN to costs corresponding to distance > searchRadius
costMat(costMat>searchRadius) = NaN;

%append matrix to allow birth and death
maxCost = max(max(max(costMat))+1,1);
deathCost = maxCost*ones(numFeaturesFrame1,1);
birthCost = maxCost*ones(numFeaturesFrame2,1);

%generate upper right and lower left block
deathBlock = diag(deathCost); %upper right
deathBlock(deathBlock==0) = NaN;
birthBlock = diag(birthCost); %lower left
birthBlock(birthBlock==0) = NaN;

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);
lrBlock = costMat';
lrBlock(~isnan(lrBlock)) = costLR;

%append cost matrix
costMat = [costMat deathBlock; birthBlock lrBlock];

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;


%%%%% ~~ the end ~~ %%%%%

%% old snippets of code

%             .cost2probExpo      : Exponent coefficient in formula
%                                   converting linking cost to scheme
%                                   probability.
% cost2probExpo    = costMatParam.cost2probExpo;
% %reserve memory
% costMat = NaN*ones(numFeaturesFrame1,numFeaturesFrame2);
% propagationScheme = ones(numFeaturesFrame1,numFeaturesFrame2);
% 
% %pick a set of uniformly distributed random numbers
% randomNumber = rand(numFeaturesFrame1,numFeaturesFrame2);
% 
% %go over all features in frame 1
% for i1 = 1 : numFeaturesFrame1
% 
%     %get the square search radius of feature i1 in frame 1
%     %     searchRadius2 = min(noiseStdMult2*2*...
%     %         kalmanFilterInfoFrame1.noiseVar(1,1,i1),maxSearchRadius2);
%     searchRadius2 = min(noiseStdMult2*2*...
%         max(kalmanFilterInfoFrame1.noiseVar(1,1,:)),maxSearchRadius2);
% 
%     %go over all features in frame 2
%     for i2 = 1 : numFeaturesFrame2
% 
%         %get the costs for this pair
%         costPair = squeeze(costMatTmp(i1,i2,:));
% 
%         %assign NaN to schemes whose cost = square distance > searchRadius2
%         costPair(costPair>searchRadius2) = NaN;
% 
%         %find the propagation scheme with the minimum cost
%         [costTmp,propSchemeTmp] = min(costPair);
% 
%         %if there is a possible link between the two features
%         if ~isnan(costTmp)
% 
%             %if Brownian motion has the minumum cost, accept it
%             if propSchemeTmp == 3
% 
%                 propagationScheme(i1,i2) = propSchemeTmp;
%                 costMat(i1,i2) = costTmp;
% 
%             else %if forward or backward propagation has minimum cost
% 
%                 %make the other propagation scheme impossible
%                 if propSchemeTmp == 1
%                     costPair(2) = NaN;
%                 else
%                     costPair(1) = NaN;
%                 end
% 
%                 %calculate the probability of each scheme using its cost
%                 schemeProb = exp(-cost2probExpo*costPair);
% 
%                 %assign a probability of zero to impossible schemes
%                 schemeProb(isnan(schemeProb)) = 0;
% 
%                 %normalize the probabilities
%                 schemeProb = schemeProb/sum(schemeProb);
% 
%                 %calculate the cumulative probability
%                 schemeProb(2) = schemeProb(2) + schemeProb(1);
%                 schemeProb(3) = 1;
% 
%                 %choose propagation scheme based on probabilities and random number
%                 propagationScheme(i1,i2) = find((randomNumber(i1,i2)<=schemeProb),1,'first');
% 
%                 %assign cost based on chosen propagation scheme
%                 costMat(i1,i2) = costMatTmp(i1,i2,propagationScheme(i1,i2));
% 
%             end %(if propSchemeTmp == 3)
% 
%         end %(if any(~isnan(costPair)))
% 
%     end %(for i2 = 1 : numFeaturesFrame2)
% end %(for i1 = 1 : numFeaturesFrame1)


% %replicate the feature amplitudes at the 2 time points to get n-by-m
% %matrices
% a1 = repmat(movieInfo(1).amp(:,1),1,numFeaturesFrame2);
% a2 = repmat(movieInfo(2).amp(:,1)',numFeaturesFrame1,1);
%
% %divide the larger of the two amplitudes by the smaller value
% ampRatio = a1./a2;
% for j=1:numFeaturesFrame2
%     for i=1:numFeaturesFrame1
%         if ampRatio(i,j) < 1
%             ampRatio(i,j) = 1/ampRatio(i,j);
%         end
%     end
% end
%
% %assign NaN to all pairs whose amplitude ratio is larger than the
% %maximum allowed
% ampRatio(ampRatio>maxAmpRatio) = NaN;
%
% %multiply the distance between pairs with the ratio between their
% %amplitudes
% costMat = costMat.*ampRatio;


% % %for features in first frame that are linked to features in previous
% % %frames, assign death cost as the 80th percentile of the cost of previous
% % %links
% % deathCost = prctile(prevCost(1:numFeaturesFrame1,:),50,2);
% % 
% % %calculate 80th percentile of all previous costs to assign as a death cost
% % %for features in 1st frame not linked to the past and as a birth cost for
% % %features in 2nd frame
% % prevCostAve = prctile(prevCost(:),50);
% % prevCostAve(isnan(prevCostAve)) = max(max(costMat)) + 1;
% % deathCost(isnan(deathCost)) = prevCostAve;
% % birthCost = prevCostAve*ones(numFeaturesFrame2,1);


