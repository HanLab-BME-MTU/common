function [trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
    linkFeaturesTp2Tp(movieInfo,searchRadius,maxAmpRatio)
%LINKFEATURESTP2TP links features between consecutive time points in a movie using LAP
%
%SYNOPSIS function [trackedFeatureNum,trackedFeatureInfo,errFlag] = ...
%    linkFeaturesTp2Tp(movieInfo,searchRadius)
%
%INPUT  movieInfo   : Array of size equal to the number of time points
%                     in a movie, containing the fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coorsinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       searchRadius: Maximum distance between two features in two
%                     consecutive time points that allows linking them (in pixels).
%       maxAmpRatio : Maximum ratio between the amplitudes of two features
%                     in two censecutive time points that allows linking them.
%
%OUTPUT trackedFeatureNum: Connectivity matrix of features between time points.
%                          Rows indicate continuous tracks, while columns 
%                          indicate time points. A track that ends before the
%                          last time point is followed by zeros, and a track
%                          that starts at a time after the first time point
%                          is preceded by zeros. 
%       trackedFeatureInfo:The positions and amplitudes of the tracked
%                          features. Number of rows = number of tracks, 
%                          while number of columns = 6*number of time 
%                          points. Each row consists of 
%                          [x1 y1 a1 dx1 dy1 da1 x2 y2 a2 dx2 dy2 da2 ...]
%                          in image coordinate system (coordinates in
%                          pixels). NaN is used to indicate time points 
%                          where the track does not exist.
%       errFlag          : 0 if function executes normally, 1 otherwise.
%
%REMARKS No gap closing in this code.
%        The cost of a link between 2 features in 2 consecutive frames
%        = (displacement between frames) times (ratio of larger amplitude
%        to smaller amplitude).
%        The algorithm is currently for the special case of 2D, but in
%        principle it can be generalized to ND quite easily.
%
%Khuloud Jaqaman, August 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackedFeatureNum = [];
trackedFeatureInfo = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('linkFeaturesTp2Tp')
    disp('--linkFeaturesTp2Tp: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of time points in movie
numTimePoints = length(movieInfo);

%get number of features at each time point
for t=1:numTimePoints
    movieInfo(t).num = size(movieInfo(t).xCoord,1);
end

%fill the feature numbers in first time point in the connectivity matrix
trackedFeatureNum = [1:movieInfo(1).num]';

%go over all time points
for t = 1:numTimePoints-1

    %calculate cost matrix, where costMat(i,j) = (distance between feature i
    %in time point t and feature j in time point t+1) times (amplitude ratio, 
    %where the larger amplitude is divided by the smaller amplitude).

    %get number of features in the 2 time points
    n = movieInfo(t).num;
    m = movieInfo(t+1).num;

    %replicate x,y-coordinates at the 2 time points to get n-by-m matrices
    x1 = repmat(movieInfo(t).xCoord(:,1),1,m);
    y1 = repmat(movieInfo(t).yCoord(:,1),1,m);
    x2 = repmat(movieInfo(t+1).xCoord(:,1)',n,1);
    y2 = repmat(movieInfo(t+1).yCoord(:,1)',n,1);

    %calculate the square distances between features in time points t and t+1
    costMat = (x1-x2).^2 + (y1-y2).^2;

    %assign NaN to all pairs that are separated by a distance > searchRadius
    indx = find(costMat>searchRadius^2);
    costMat(indx) = NaN;

    %replicate the feature amplitudes at the 2 time points to get n-by-m
    %matrices
    a1 = repmat(movieInfo(t).amp(:,1),1,m);
    a2 = repmat(movieInfo(t+1).amp(:,1)',n,1);
    
    %divide the larger of the two amplitudes by the smaller value
    ampRatio = a1./a2;
    for j=1:m
        for i=1:n
            if ampRatio(i,j) < 1
                ampRatio(i,j) = 1/ampRatio(i,j);
            end
        end
    end    
    
    %assign NaN to all pairs whose amplitude ratio is larger than the
    %maximum allowed
    indx = find(ampRatio>maxAmpRatio);
    ampRatio(indx) = NaN;
    
    clear x1 y1 x2 y2 a1 a2;

    %multiply the distance between pairs with the ratio between their
    %amplitudes
    costMat = costMat.*ampRatio;

    %replace NaN, indicating pairs separated by a distance > searchRadius,
    %with -1
    indx = find(isnan(costMat));
    costMat(indx) = -1;

    %track features based on this cost matrix, allowing for birth and death
    [link12,link21] = lap(costMat,-1,0,1);
    
    %get indices of features at time t+1 that are connected to features at time t
    indx2C = find(link21(1:m)<=n);
    
    %get indices of corresponding features at time t
    indx1C = link21(indx2C);
    
    %find the rows in "trackedFeatureNum" that are not connected to features at time t+1
    indx1U = ones(size(trackedFeatureNum,1),1);
    indx1U(indx1C) = 0;
    indx1U = find(indx1U);
    
    %assign space for new matrix
    tmp = zeros(size(trackedFeatureNum,1)+m-length(indx2C),t+1);
    
    %fill in the feature numbers at time t+1
    tmp(1:m,t+1) = [1:m]';
    
    %shuffle the rows from the previous times to get the correct
    %connectivity with time point t+1
    tmp(indx2C,1:t) = trackedFeatureNum(indx1C,:);
    
    %add rows of tracks that are not connected to points at time t+1
    tmp(max(m,length(indx1C))+1:end,1:t) = trackedFeatureNum(indx1U,:);

    %update the connectivity matrix "trackedFeatureNum"
    trackedFeatureNum = tmp;
    
end

%get total number of tracks
numTracks = size(trackedFeatureNum,1);

%find the time point where each track begins and then sort the vector
tpStart = zeros(numTracks,1);
for i=1:numTracks
    tpStart(i) = find((trackedFeatureNum(i,:)~=0),1,'first');
end
[tpStart,indx] = sort(tpStart);

%rearrange "trackedFeatureNum" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the 
%same time in descending order from longest to shortest.
trackedFeatureNum = trackedFeatureNum(indx,:);

%store feature positions and amplitudes in a matrix that also shows their connectivities
%information is stored as as [x y a dx dy da] in image coordinate system
numRows = size(trackedFeatureNum,1);
trackedFeatureInfo = NaN*ones(numRows,6*numTimePoints);
for t=1:numTimePoints
    indx1 = find(trackedFeatureNum(:,t)~=0);
    indx2 = trackedFeatureNum(indx1,t);
    trackedFeatureInfo(indx1,6*(t-1)+1:6*t) = [movieInfo(t).xCoord(indx2,1) ...
        movieInfo(t).yCoord(indx2,1) movieInfo(t).amp(indx2,1) ...
        movieInfo(t).xCoord(indx2,2) movieInfo(t).yCoord(indx2,2) ...
        movieInfo(t).amp(indx2,2)];
end


%%%%% ~~ the end ~~ %%%%%
