function [trajectoryOut, linFit] = removeLinearTrend(trajectory,sigmaTrajectory,robust)
%REMOVELINEARTREND removes a linear trend from a data series 
%   by doing a robust linear fit and subsequent subtraction. Furthermore,
%   it sets the mean to zero.
%
% SYNOPSIS  outTrajectory = removeLinearTrend(trajectory, sigmaTrajectory, robust)
%
% INPUT     trajectory     : any 1D data series vector. Can contain NaNs
%           sigmaTrajectory: (optional) The corresponding uncertainties
%           robust         : (optional) If true, a least-median squares fit
%                             will be used to de-trend. If false, a
%                             least-squares fit will be used. Default is
%                             true.
%
% OUTPUT    trajectoryOut: the above data series without linear trend and
%                          with zero mean
%           linFit       : intercept and slope of linear fit
%
%
% currently, there is no error propagation implemented.
%
%c: jonas 9/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%============
% TEST INPUT
%============

% get size of original trajectory
%trajSizeOri = size(trajectory);

% make into column vector
trajectory = returnRightVector(trajectory);

% and get the new size
trajSize = size(trajectory);

% check for sigma
if nargin > 1 && ~isempty(sigmaTrajectory)
    sigmaTrajectory = returnRightVector(sigmaTrajectory);
else
    sigmaTrajectory = ones(trajSize);
end

%Check for robust
if nargin < 3 || isempty(robust)
    robust = true;
end
    
%============


%=================
% LINEAR FIT
%=================

if robust
    % we fit A*x=B+E

    % prepare matrices
    A = [ones(trajSize),[1:trajSize(1)]'];
    B = trajectory;

    % robust linear fit
    linFit = linearLeastMedianSquares(A,B,diag(sigmaTrajectory));
    
    %=================

    %===================
    % SUBTRACT & RETURN
    %===================

    % subtract trend
    trajectory = trajectory - A*linFit;

else%Added option to do non-robust fit (it's MUCH faster) - Hunter 3/2011 
    
    %Deal with NaNs
    x = (1:trajSize)';
    y = trajectory;
    iNan = isnan(y);
    x(iNan) = [];
    y(iNan) = [];
    %Do least-squares fit
    p = polyfit(x,y,1);
    %Subtract trend
    trajectory = trajectory - polyval(p,1:trajSize)';        
    
    
    
end
    


% subtract mean
trajectoryOut = trajectory - nanmean(trajectory);
