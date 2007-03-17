function [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo)
%KALMANINITLINEARMOTION initializes Kalman filter state vector and covariance matrix for features in one frame
%
%SYNOPSIS [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo)
%
%INPUT  frameInfo       : Structure with fields (for 1 frame):
%             .xCoord       : Image coordinate system x-coordinates of detected
%                             features (in pixels). 1st column for
%                             value and 2nd column for standard deviation.
%             .yCoord       : Image coordinate system y-coordinates of detected
%                             features (in pixels). 1st column for
%                             value and 2nd column for standard deviation.
%             .amp          : Amplitudes of PSFs fitting detected features.
%                             1st column for values and 2nd column
%                             for standard deviations.
%             .num          : Number of features in frame.
%
%OUTPUT kalmanFilterInfo: Structure with fields:
%             .stateVec     : State vector for each feature.
%             .stateCov     : State covariance matrix for each feature.
%             .noiseVar     : Variance of state noise for each feature.
%       errFlag         : 0 if function executes normally, 1 otherwise.
%
%
%Khuloud Jaqaman, March 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kalmanFilterInfo = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--kalmanInitLinearMotion: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find number of features in frame
numFeatures = frameInfo.num;

%initialize state vector
kalmanFilterInfo.stateVec = [frameInfo.xCoord(:,1) ...
    frameInfo.yCoord(:,1) zeros(numFeatures,2)];

%initialize state covariance matrix
for iFeature = 1 : numFeatures
    kalmanFilterInfo.stateCov(:,:,iFeature) = diag(...
        [frameInfo.xCoord(iFeature,2)^2 frameInfo.yCoord(iFeature,2)^2 4 4]);
end

%initialize state noise variance
for iFeature = 1 : numFeatures
    kalmanFilterInfo.noiseVar(:,:,iFeature) = diag(...
        [4 4 4 4]);
end


%%%%% ~~ the end ~~ %%%%%
