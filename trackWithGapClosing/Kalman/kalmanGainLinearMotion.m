function [kalmanFilterInfoOut,errFlag] = kalmanGainLinearMotion(trackedFeatureIndx,...
    frameInfo,kalmanFilterInfoTmp,propagationScheme,kalmanFilterInfoIn,...
    filterInfoPrev,kalmanInitParam)
%KALMANGAINLINEARMOTION uses the Kalman gain from linking to get better estimates of the state vector, its covariance matrix, state noise and its variance
%
%SYNOPSIS [kalmanFilterInfoOut,errFlag] = kalmanGainLinearMotion(trackedFeatureIndx,...
%    frameInfo,kalmanFilterInfoTmp,propagationScheme,kalmanFilterInfoIn,...
%    filterInfoPrev,kalmanInitParam)
%
%INPUT  trackedFeatureIndx : Matrix showing connectivity between features
%                            in current frame (listed in last column of matrix) 
%                            and features in previous frames. A zero in 
%                            columns before last indicates that feature is
%                            not connected to any previous features.
%       frameInfo          : Structure with fields (for current frame):
%             .xCoord          : Image coordinate system x-coordinates of detected
%                                features (in pixels). 1st column for
%                                value and 2nd column for standard deviation.
%             .yCoord          : Image coordinate system y-coordinates of detected
%                                features (in pixels). 1st column for
%                                value and 2nd column for standard deviation.
%             .amp             : Amplitudes of PSFs fitting detected features. 
%                                1st column for values and 2nd column 
%                                for standard deviations.
%       kalmanFilterInfoTmp: Structure with fields (for current frame):
%             .stateVec        : Kalman filter prediction of state
%                                vector of all features in current frame 
%                                based on all 3 motion models.
%             .stateCov        : Kalman filter prediction of state
%                                covariance matrix of all features in
%                                current frame based on all 3 motion models.
%             .obsVec          : Kalman filter prediction of the
%                                observed variables for all features in 
%                                current frame based on all 3 motion models.
%       propagationScheme  : Matrix indicating the propagation scheme that
%                            yielded the lowest cost for a link between two
%                            features.
%       kalmanFilterInfoIn : Structure with fields (for all previous frames):
%             .stateVec        : Kalman filter state vector for all features.
%             .stateCov        : Kalman filter state covariance matrix
%                                for all features.
%             .noiseVar        : Variance of state noise for all feature.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       filterInfoPrev     : Structure with fields (for current frame):
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%                            Optional input. Enter [] or nothing if not to be used.
%       kalmanInitParam    : Structure with fields containing variables
%                            used in Kalman filter initialization. See
%                            particular initialization function for fields.
%                            Optional. Enter [] or nothing if not to be used.
%
%OUTPUT kalmanFilterInfoOut: Structure with fields (for all frames upto current):
%             .stateVec        : Kalman filter state vector for all features.
%             .stateCov        : Kalman filter state covariance matrix
%                                for all features.
%             .noiseVar        : Variance of state noise for all features.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       errFlag            : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kalmanFilterInfoOut = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 5
    disp('--kalmanGainLinearMotion: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether a priori Kalman filter information is given
if nargin < 6 || isempty(filterInfoPrev)
    filterInfoPrev = [];
    usePriorInfo = 0;
else
    usePriorInfo = 1;
end

%check whether additional parameters for Kalman filter initialization are
%given
if nargin < 7 || isempty(kalmanInitParam)
    kalmanInitParam = [];
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gain calculation and update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%copy kalmanFilterInfoIn into kalmanFilterInfoOut
kalmanFilterInfoOut = kalmanFilterInfoIn;

%get number of features in current frame and frame number
[numFeatures,iFrame] = size(trackedFeatureIndx);

%construct Kalman filter observation matrix
observationMat = [1 0 0 0; 0 1 0 0];

%go over all features in current frame
for iFeature = 1 : numFeatures

    %find index of feature in previous frame that this feature is connected to
    iFeaturePrev = trackedFeatureIndx(iFeature,end-1);

    %if this feature is connected to a feature in previous frame
    if iFeaturePrev ~= 0

        %find propagation scheme leading to this link and save in kalmanFilterInfo
        iScheme = propagationScheme(iFeaturePrev,iFeature);
        kalmanFilterInfoOut(iFrame).scheme(iFeature,1) = iScheme; %to current feature
        kalmanFilterInfoOut(iFrame-1).scheme(iFeaturePrev,2) = iScheme; %from previous feature

        %get the corresponding state vector, covariance and observation vector
        stateVecOld = kalmanFilterInfoTmp.stateVec(iFeaturePrev,:,iScheme)';
        stateCovOld = kalmanFilterInfoTmp.stateCov(:,:,iFeaturePrev,iScheme);
        obsVecOld = kalmanFilterInfoTmp.obsVec(iFeaturePrev,:,iScheme)';
        
        %calculate Kalman gain
        kalmanGain = stateCovOld*observationMat'*...
            inv(observationMat*stateCovOld*observationMat'+...
            diag([frameInfo.xCoord(iFeature,2)^2 frameInfo.yCoord(iFeature,2)^2]));

        %estimate state noise in previous frame and save in kalmanFilterInfo
        stateNoise = kalmanGain*([frameInfo.xCoord(iFeature,1) ...
            frameInfo.yCoord(iFeature,1)]'-obsVecOld);
        kalmanFilterInfoOut(iFrame-1).stateNoise(iFeaturePrev,:) = stateNoise';
        
        %update estimate of state vector in current frame
        stateVec = stateVecOld + stateNoise;
        
        %update estimate of state covariance matrix in current frame
        stateCov = stateCovOld - kalmanGain*observationMat*stateCovOld;

        %find indices of previous features connected to this feature
        indx = trackedFeatureIndx(iFeature,1:end-1); 
        
        %determine where the track starts
        indxLength = length(find(indx));
        
        %collect all of the error terms
        stateNoiseAll = [];
        for i = iFrame-indxLength : iFrame-1
            stateNoiseAll = [stateNoiseAll; kalmanFilterInfoOut(i).stateNoise(indx(i),:)];
        end

        %impose isotropy (x and y are equivalent)
        stateNoisePos = [stateNoiseAll(:,1); stateNoiseAll(:,2)];
        stateNoiseVel = [stateNoiseAll(:,3); stateNoiseAll(:,4)];

        %         %don't impose isotropy
        %         stateNoisePos = stateNoiseAll(:,1:2);
        %         stateNoiseVel = stateNoiseAll(:,3:4);

        %estimate positional noise variance in current frame
        noiseVar(1:2) = var(stateNoisePos);

        %estimate speed noise variances in current frame
        noiseVar(3:4) = var(stateNoiseVel);

        %save this information in kalmanFilterInfo
        kalmanFilterInfoOut(iFrame).stateVec(iFeature,:) = stateVec';
        kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature) = stateCov;
        kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature) = diag(noiseVar);

    else %if this feature is not connected to anything in previous frame

        %initialize Kalman filter for this feature
        if usePriorInfo %use prior information if supplied
            kalmanFilterInfoOut(iFrame).stateVec(iFeature,:) = filterInfoPrev.stateVec(iFeature,:);
            kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature) = filterInfoPrev.stateCov(:,:,iFeature);
            kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature) = filterInfoPrev.noiseVar(:,:,iFeature);
        else
            featureInfo.xCoord = frameInfo.xCoord(iFeature,:);
            featureInfo.yCoord = frameInfo.yCoord(iFeature,:);
            featureInfo.amp = frameInfo.amp(iFeature,:);
            featureInfo.num = 1;
            [filterTmp,errFlag] = kalmanInitLinearMotion(featureInfo,kalmanInitParam);
            kalmanFilterInfoOut(iFrame).stateVec(iFeature,:) = filterTmp.stateVec;
            kalmanFilterInfoOut(iFrame).stateCov(:,:,iFeature) = filterTmp.stateCov;
            kalmanFilterInfoOut(iFrame).noiseVar(:,:,iFeature) = filterTmp.noiseVar;
        end

    end %(if iFeaturePrev ~= 0 ... else ...)

end %(for iFeature = 1 : numFeatures)


%%%%% ~~ the end ~~ %%%%%
