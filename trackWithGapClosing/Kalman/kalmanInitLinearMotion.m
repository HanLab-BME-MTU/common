function [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo,...
    probDim,initParam)
%KALMANINITLINEARMOTION initializes Kalman filter state vector and covariance matrix for features in one frame
%
%SYNOPSIS [kalmanFilterInfo,errFlag] = kalmanInitLinearMotion(frameInfo,...
%    probDim,initParam)
%
%INPUT  frameInfo       : Structure with fields (for 1 frame):
%             .allCoord     : Image coordinate system x,dx,y,dy,[z,dz] of 
%                             detected features and their uncertainties (in
%                             pixels).
%             .num          : Number of features in frame.
%       probDim         : Problem dimension. 2 (for 2D) or 3 (for 3D).
%       initParam       : Structure with fields
%             .convergePoint: Convergence point (x, y, [z] coordinates) of tracks
%                             if motion is radial, in image coordinate
%                             system. Should be a row vector.
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
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--kalmanInitLinearMotion: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check whether convergence point for radial motion is supplied
if nargin < 3 || isempty(initParam) || ~isfield(initParam,'convergePoint')
    initParam.convergePoint = [];
end

%get convergence point for radial motion
convergePoint = initParam.convergePoint;

if ~isempty(convergePoint)
    [nRow,nCol] = size(convergePoint);
    if nRow ~= 1 && nCol == 1
        convergePoint = convergePoint';
    end
    nCol = size(convergePoint,2);
    if nCol == probDim
        disp('--kalmanInitLinearMotion: initParam.convergePoint of wrong dimension!');
        errFlag = 1;
    end
end
    
%exit if there are problems with input
if errFlag
    disp('--kalmanInitLinearMotion: Please fix input parameters.');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find number of features in frame
numFeatures = frameInfo.num;

%estimate initial velocity
if isempty(convergePoint) %if there is no convergence point information

    %assume zero initial velocity
    velocityInit = zeros(numFeatures,probDim);

else %if there is a convergence point

    %assign initial speed
    speedInit = 1;

    %get displacement and distance of features from convergence point
    displacement = repmat(convergePoint,numFeatures,1) - ...
        frameInfo.allCoord(1:2:end);
    distance = sqrt(sum(displacement.^2,2));

    %calculate initial velocity
    velocityInit = speedInit * displacement ./ repmat(distance,1,probDim);

end

%initialize state vector
kalmanFilterInfo.stateVec = [frameInfo.allCoord(:,1:2:end) velocityInit];

%initialize state covariance matrix
for iFeature = 1 : numFeatures
    kalmanFilterInfo.stateCov(:,:,iFeature) = diag(...
        [frameInfo.allCoord(iFeature,2:2:end).^2 4*ones(1,probDim)]);
end

%initialize state noise variance
for iFeature = 1 : numFeatures
    kalmanFilterInfo.noiseVar(:,:,iFeature) = diag(...
        ones(2*probDim,1));
end


%%%%% ~~ the end ~~ %%%%%
