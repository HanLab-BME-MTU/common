function [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker,...
     errFlag] = costMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,costMatParam)
%COSTMATLINEARMOTIONLINK provides a cost matrix for linking features based on competing linear motion models
%
%SYNOPSIS [costMat,propagationScheme,kalmanFilterInfoFrame2,nonlinkMarker] ...
%     = costMatLinearMotionLink(movieInfo,kalmanFilterInfoFrame1,costMatParam)
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
%      kalmanFilterInfoFrame1 : Structure with the following fields:
%             .stateVec           : Kalman filter state vector for each
%                                   feature in 1st frame.
%             .stateCov           : Kalman filter state covariance matrix
%                                   for each feature in 1st frame.
%             .noiseVar           : Variance of state noise for each
%                                   feature in 1st frame.
%      costMatParam           : Structure with fields:
%             .maxSearchRadius    : Maximum allowed search radius (in pixels).
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
%        1: zero drift (Brownian), 2: forward drift, 3: backward drift.
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
if nargin ~= nargin('costMatLinearMotion')
    disp('--costMatLinearMotion: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

maxSearchRadius2 = (costMatParam.maxSearchRadius)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify number of propagation schemes used
numSchemes = 3;

%construct transition matrices
transMat(:,:,1) = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1]; %forward drift transition matrix
transMat(:,:,2) = [1 0 -1 0; 0 1 0 -1; 0 0 1 0; 0 0 0 1]; %backward drift transition matrix
% transMat(:,:,1) = eye(4); %zero drift transition matrix
% transMat(:,:,2) = eye(4); %zero drift transition matrix
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
for iFeature = 1:numFeaturesFrame1

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

%replicate the observed positions in the 2nd frame to get 
%numFeaturesFrame1 x numFeaturesFrame2 matrices
x2 = repmat(movieInfo(2).xCoord(:,1)',numFeaturesFrame1,1);
y2 = repmat(movieInfo(2).yCoord(:,1)',numFeaturesFrame1,1);

%calculate the cost matrices for all three propagation scheme
for iScheme = 1 : numSchemes

    %replicate the forward propagated positions in the 1st frame
    x1 = repmat(propagatedPos(:,1,iScheme),1,numFeaturesFrame2);
    y1 = repmat(propagatedPos(:,2,iScheme),1,numFeaturesFrame2);

    %calculate the square distances between features
    costMat(:,:,iScheme) = (x1-x2).^2 + (y1-y2).^2;

end

%find the minimum cost for the link between every pair, which also 
%determines the best propagation scheme to perform that link
[costMat,propagationScheme] = min(costMat,[],3);

%estimate square search radius from positional noise variance
searchRadius2 = min(25*max(kalmanFilterInfoFrame1.noiseVar(1,1,:)),maxSearchRadius2);

%assign NaN to all pairs that are separated by a distance > searchRadius
costMat(costMat>searchRadius2) = NaN;

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

%determine the nonlinkMarker
nonlinkMarker = min(floor(min(min(costMat)))-5,-5);

%replace NaN, indicating pairs that cannot be linked, with nonlinkMarker
costMat(isnan(costMat)) = nonlinkMarker;


%%%%% ~~ the end ~~ %%%%%

