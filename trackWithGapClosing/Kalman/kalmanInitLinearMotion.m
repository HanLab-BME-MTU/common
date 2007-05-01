function [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo,initParam)
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
%       initParam       : Structure with fields
%             .convergePoint: Convergence point (x and y coordinates) of tracks
%                             if motion is radial, in image coordinate system.
%                             Optional. If supplied, radial form is assumed and
%                             value is used to estimate initial velocities. If
%                             not supplied, then initial velocity is taken as
%                             zero. Default: [].
%                         Optional. Enter as [] if not used.
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

%check whether convergence point for radial motion is supplied
if nargin < 2 || isempty(initParam) || ~isfield(initParam,'convergePoint')
    initParam.convergePoint = [];
end

%get convergence point for radial motion
convergePoint = initParam.convergePoint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find number of features in frame
numFeatures = frameInfo.num;

%estimate initial velocity
if isempty(convergePoint) %if there is no convergence point information

    %assume zero initial velocity
    velocityInit = zeros(numFeatures,2);

else %if there is a convergence point

    %assign initial speed
    speedInit = 1;

    %get displacements in x and y and distance between features and convergence point
    xDisp = convergePoint(1) - frameInfo.xCoord(:,1);
    yDisp = convergePoint(2) - frameInfo.yCoord(:,1);
    distance = sqrt(sum([xDisp yDisp].^2,2));

    %calculate initial velocity
    velocityInit = speedInit*[xDisp./distance yDisp./distance];

end

%initialize state vector
kalmanFilterInfo.stateVec = [frameInfo.xCoord(:,1) ...
    frameInfo.yCoord(:,1) velocityInit];

%initialize state covariance matrix
for iFeature = 1 : numFeatures
    kalmanFilterInfo.stateCov(:,:,iFeature) = diag(...
        [frameInfo.xCoord(iFeature,2)^2 frameInfo.yCoord(iFeature,2)^2 4 4]);
end

%initialize state noise variance
for iFeature = 1 : numFeatures
    kalmanFilterInfo.noiseVar(:,:,iFeature) = diag(...
        [1 1 1 1]);
end


%%%%% ~~ the end ~~ %%%%%
