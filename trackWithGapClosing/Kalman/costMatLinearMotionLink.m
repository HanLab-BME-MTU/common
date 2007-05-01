function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
     errFlag] = costMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,...
     costMatParam,useLocalDensity,nnDistTracks)
%COSTMATLINEARMOTIONLINK provides a cost matrix for linking features based on competing linear motion models
%
%SYNOPSIS [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
%     errFlag] = costMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,...
%     costMatParam,useLocalDensity,nnDistTracks)
%
%INPUT  movieInfo             : A 2x1 array (corresponding to the 2 frames of 
%                               interest) containing the fields:
%             .xCoord             : Image coordinate system x-coordinates of detected
%                                   features (in pixels). 1st column for
%                                   value and 2nd column for standard deviation.
%             .yCoord             : Image coordinate system y-coordinates of detected
%                                   features (in pixels). 1st column for
%                                   value and 2nd column for standard deviation.
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
%             .maxSearchRadius    : Maximum allowed search radius (in pixels).
%             .minSearchRadius    : Minimum allowed search radius (in pixels).
%             .brownStdMultLink   : Factor multiplying Brownian
%                                   displacement std to get search radius.
%             .closestDistScaleLink:Scaling factor of nearest neighbor
%                                   distance. Not needed if useLocalDensity = 0;
%             .maxStdMultLink     : Maximum value of factor multiplying
%                                   Brownian displacement std to get search
%                                   radius. Not needed if useLocalDensity = 0;
%      useLocalDensity        : Logical variable indicating whether to use
%                               local density in search radius estimation.
%      nnDistTracks           : Nearest neighbor distance of features in
%                               frame 1 given their history.
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
maxSearchRadius = costMatParam.maxSearchRadius;
minSearchRadius = costMatParam.minSearchRadius;
brownStdMult    = costMatParam.brownStdMultLink;
if useLocalDensity
    closestDistScale = costMatParam.closestDistScaleLink;
    maxStdMult = costMatParam.maxStdMultLink;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify number of propagation schemes used
numSchemes = 3;

%construct transition matrices
transMat(:,:,1) = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; %forward drift transition matrix
transMat(:,:,2) = [1 0 -1 0; 0 1 0 -1; 0 0 1 0; 0 0 0 1]; %backward drift transition matrix
transMat(:,:,3) = eye(4); %zero drift transition matrix

%construct observation matrix
observationMat = [1 0 0 0; 0 1 0 0]; %observation matrix

%get number of features in the 2 frames
numFeaturesFrame1 = movieInfo(1).num;
numFeaturesFrame2 = movieInfo(2).num;

%reserve memory for "kalmanFilterInfoframe2"
kalmanFilterInfoFrame2 = struct('stateVec',zeros(numFeaturesFrame1,4,numSchemes),...
    'stateCov',zeros(4,4,numFeaturesFrame1,numSchemes),...
    'obsVec',zeros(numFeaturesFrame1,2,numSchemes));

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

%put the x and y coordinates of features in the 2nd frame in one matrix
yx2 = [movieInfo(2).yCoord(:,1) movieInfo(2).xCoord(:,1)];

%calculate the cost matrices for all three propagation schemes
for iScheme = 1 : numSchemes

    %put the propagated x and y coordinates of features from 1st frame in
    %one matrix
    yx1 = [propagatedPos(:,2,iScheme) propagatedPos(:,1,iScheme)];

    %calculate the distances between features
    costMatTmp(:,:,iScheme) = createDistanceMatrix(yx1,yx2);

end

%find the minimum cost for the link between every pair, which also
%determines the best propagation scheme to perform that link
[costMat,propagationScheme] = min(costMatTmp,[],3);

%get the Kalman standard deviation of all features in frame 1
kalmanStd = sqrt(2*squeeze(kalmanFilterInfoFrame1.noiseVar(1,1,:)));

%copy brownStdMult into vector
stdMultInd = repmat(brownStdMult,numFeaturesFrame1,1);

%if local density information is used to expand search radius ...
if useLocalDensity

    %divide each feature's nearest neighbor distance/closestDistScale by kalmanStd
    ratioDist2Std = nnDistTracks./kalmanStd/closestDistScale;
   
    %make ratios larger than maxStdMult equal to 10
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

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;


%%%%% ~~ the end ~~ %%%%%


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



