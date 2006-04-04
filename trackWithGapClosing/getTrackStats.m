function [trackStats,statsRelChange,errFlag] = getTrackStats(trackedFeatureInfo,lenFrac,...
    timeWindow,trackStatsOld)
%GETTRACKSTATS determines the statistical characeteristics of amplitude change and displacement in tracks over time
%
%SYNOPSIS [trackStats,errFlag] = getTrackStats(trackedFeatureInfo,lenFrac,...
%    timeWindow)
%
%INPUT  trackedFeatureInfo:The positions and amplitudes of the tracked
%                          features. Number of rows = number of tracks, 
%                          while number of columns = 6*number of time 
%                          points. Each row consists of 
%                          [x1 y1 a1 dx1 dy1 da1 x2 y2 a2 dx2 dy2 da2 ...]
%                          in image coordinate system (coordinates in
%                          pixels). NaN is used to indicate time points 
%                          where the track does not exist.
%       lenFrac          : Minimum length of tracks used for statistical
%                          analysis, as a fraction of the total number of
%                          time points in movie.
%       timeWindow       : Time window of gap closing.
%                          Optional. Default: 1.
%       trackStatsOld    : trackStats (See output description) from previous 
%                          calculation. Optional. Default: [].
%
%OUTPUT trackStats       : Structure with fields:
%           .dispSqLambda     : timeWindow x 1 vector of parameter of the 
%                               exponential distribution that describes the 
%                               displacement of a feature between frames.
%           .ampDiffStd       : timeWindow x 1 vector of standard deviation of 
%                               the change in a feature's amplitude between
%                               frames.
%       statsRelChange   : Relative change in statistical parameters.
%       errFlag          : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackStats = [];
statsRelChange = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--getTrackStats: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%assign defaults
timeWindow_def = 1;

%check timeWindow
if nargin < 3 || isempty(timeWindow)
    timeWindow = timeWindow_def;
end

if nargin < 4 || isempty(trackStatsOld)
    trackStatsOld = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of time points and tracks in movie
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints = numTimePoints/6;

%find tracks that contain more than (lenFrac*numTimePoints) time points
trackLength = zeros(numTracks,1);
for i=1:numTracks
    trackLength(i) = length(find(~isnan(trackedFeatureInfo(i,:))))/6;
end
goodTracks = find(trackLength>=lenFrac*numTimePoints);

if ~isempty(goodTracks) %if there are tracks to use ...

    %get the x,y-coordinates and amplitudes over time of these good tracks
    xCoord = trackedFeatureInfo(goodTracks,1:6:end);
    yCoord = trackedFeatureInfo(goodTracks,2:6:end);
    amplitude = trackedFeatureInfo(goodTracks,3:6:end);

    %reserve memory for output vectors
    dispSqLambda = zeros(timeWindow,1);
    ampDiffStd = zeros(timeWindow,1);
    
    for i=1:timeWindow

        %get the squared displacement between time points
        dispSq = (xCoord(:,i+1:end) - xCoord(:,1:end-i)).^2 + ...
            (yCoord(:,i+1:end) - yCoord(:,1:end-i)).^2;

        %obtain the histogram of the squared displacement
        [n,x] = histogram(dispSq(:));

        %fit the histogram with an exponential function
        expParam = lsqcurvefit(@expFun,[100 1]',x,n);
        dispSqLambda(i) = expParam(2);

        %get the amplitude change between time points
        ampDiff = amplitude(:,i+1:end) - amplitude(:,1:end-i);

        %get the standard deviation of amplitude change
        ampDiffStd(i) = nanstd(ampDiff(:));

    end %(for i=1:timeWindow)
    
    %save output in structure
    trackStats.dispSqLambda = dispSqLambda;
    trackStats.ampDiffStd = ampDiffStd;
    
    %get the maximum relative change in parameters, if the old parameters
    %are supplied
    if ~isempty(trackStatsOld)
        oldParam = [trackStatsOld.dispSqLambda;trackStatsOld.ampDiffStd];
        statsRelChange = max(abs(([dispSqLambda;ampDiffStd]-oldParam)./oldParam));
    end

else %if there aren't any ...

    errFlag = 1;
    return

end %(if ~isempty(goodTracks) ... else ...)


%%%%% ~~ the end ~~ %%%%%

