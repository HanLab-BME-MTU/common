function [trackedFeatureNum,trackedFeatureInfo,errFlag] = trackWithGapClosing(...
    movieInfo,costMatrices,trackStatFun,gapCloseParam,iterParam)
%TRACKWITHGAPCLOSING links features between frames in a movie and closes gaps, with the possibility of track merging and splitting
%
%SYNOPSIS [trackedFeatureNum,trackedFeatureInfo,errFlag] = trackWithGapClosing(...
%    movieInfo,linkingLimits,gapCloseParam,iterParam)
%
%INPUT  movieInfo    : Array of size equal to the number of time points
%                      in a movie, containing the fields:
%             .xCoord      : Image coordinate system x-coordinate of detected
%                            features (in pixels).
%             .yCoord      : Image coorsinate system y-coordinate of detected
%                            features (in pixels).
%             .amp         : Amplitudes of PSFs fitting detected features.
%       costMatrices : 4-by-1 array indicating cost matrices and their
%                      parameters.
%                      -1st entry indicates the initial simple cost matrix
%                      that links between consecutive frames.
%                      -2nd entry indicates cost matrix for linking between
%                      consecutive frames based on statistical
%                      characteristics of tracked features.
%                      -3rd entry indicates cost matrix for gap closing,
%                      also based on statistical characteristics of
%                      tracked features.
%                      -4th entry indicates cost matrix for linking a track
%                      that split from another one to those that merged to
%                      the latter, based on statistical characteristics of
%                      tracked features.
%                      Each entry is a structure with the fields:
%             .costMatFun  : Name of function used to calculate cost matrix.
%             .costMatParam: Structure containing parameters needed for cost matrix.
%       trackStatFun : Name of function to be used for extracting tracked
%                      feature statistics.
%       gapCloseParam: Structure containing variables needed for gap closing.
%                      Contains the fields:
%             .timeWindow  : Largest time gap between the end of a track and the
%                            beginning of another that could be connected to it.
%             .mergeSplit  : Logical variable with value 1 if the merging
%                            and splitting of trajectories are to be consided;
%                            and 0 if merging and splitting are not allowed.
%       iterParam    : Structure with parameters related to iterating:
%             .tolerance   : Tolerance for changes in track statistics to
%                            stop iterating.
%             .lenFrac     : Minimum length of tracks used for statistical
%                            analysis, as a fraction of the total number of
%                            time points in movie.
%
%OUTPUT trackedFeatureNum : Connectivity matrix of features between time points.
%                           Rows indicate continuous tracks, while columns
%                           indicate time points. A track that ends before the
%                           last time point is followed by zeros, and a track
%                           that starts at a time after the first time point
%                           is preceded by zeros.
%       trackedFeatureInfo: The positions and amplitudes of the tracked
%                           features. Number of rows = number of tracks,
%                           while number of columns = 6*number of time
%                           points. Each row consists of
%                           [x1 y1 a1 dx1 dy1 da1 x2 y2 a2 dx2 dy2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points
%                           where the track does not exist.
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm is currently limited to the special case of 2D, but
%        in principle it can be generalized to ND quite easily.
%
%Khuloud Jaqaman, March 2006

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
if nargin ~= 5
    disp('--trackWithGapClosing: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get parameters from input
timeWindow = gapCloseParam.timeWindow;
mergeSplit = gapCloseParam.mergeSplit;
tolerance = iterParam.tolerance;
lenFrac = iterParam.lenFrac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of time points in movie
numTimePoints = length(movieInfo);

%get initial tracks by linking features between consecutive frames using
%the simple initial criteria
[trackedFeatNumInit,trackedFeatInfoInit,errFlag] = linkFeaturesTp2Tp(...
    movieInfo,costMatrices(1).costMatFun,costMatrices(1).costMatParam);

%get track statistics
eval(['[trackStatsT,statsRelChange,errFlag] = ' trackStatFun '(trackedFeatInfoInit,lenFrac,timeWindow);'])

%exit at this point if statistical analysis could not be performed
if errFlag
    disp('--trackWithGapClosing: Getting track statistics failed. Stopping after initial linking!');
    iterate = 0;
else %otherwise iterate
    iterate = 1;
end

while iterate

    %get latest statistics
    trackStats = trackStatsT;

    %link features between consecutive frames using the information obtained
    %via the statistical analysis of tracks
    costMatParams = costMatrices(2).costMatParam;
    costMatParams.trackStats = trackStats;
    [trackedFeatNumInit,trackedFeatInfoInit,errFlag] = linkFeaturesTp2Tp(...
        movieInfo,costMatrices(2).costMatFun,costMatParams);

    %close gaps using the information obtained via the statistical analysis
    %of tracks ...

    %get number of tracks formed by initial linking
    numTracks = size(trackedFeatNumInit,1);

    %find the starting time points of all tracks, find tracks that start
    %after the first time point, and get their number
    trackStartTime = zeros(numTracks,1);
    for i=1:numTracks
        trackStartTime(i) = find((trackedFeatNumInit(i,:)~=0),1,'first');
    end
    indxStart = find(trackStartTime>1);
    m = length(indxStart);

    %find the termination time points of all tracks, find tracks that end
    %before the last time point, and get the number
    trackEndTime = zeros(numTracks,1);
    for i=1:numTracks
        trackEndTime(i) = find((trackedFeatNumInit(i,:)~=0),1,'last');
    end
    indxEnd = find(trackEndTime<numTimePoints);
    n = length(indxEnd);

    %if there are gaps to close ...
    if n~=0 && m~=0
        
        %calculate the cost matrix, in sparse matrix format
        costMatParams = costMatrices(3).costMatParam;
        costMatParams.trackStats = trackStats;
        eval(['[costMat,noLinkCost,trackStartTime,trackEndTime,indxMerge,' ...
            'numMerge,indxSplit,numSplit,errFlag] =' costMatrices(3).costMatFun ...
            '(trackedFeatInfoInit,trackStartTime,indxStart,'...
            'trackEndTime,indxEnd,costMatParams,gapCloseParam);'])

        %link tracks based on this cost matrix, allowing for birth and death
        [link12,link21] = lap(costMat,-1000,0,1,noLinkCost);

        %Connect the tracks using the information from LAP
        trackedFeatureNum  = trackedFeatNumInit;
        trackedFeatureInfo = trackedFeatInfoInit;
        
        for i=m:-1:1 %go over all track starts

            if link21(i) <= n %if this start is linked to an end

                %get the track to be appended at its end, the time point at which
                %it will be appended, and the track that will be added to it
                track2Append = indxEnd(link21(i));
                time2Append = trackStartTime(i);
                trackAdded = indxStart(i);

                %modify the matrix indicating linked feature number
                trackedFeatureNum(track2Append,time2Append:end) = ...
                    trackedFeatureNum(trackAdded,time2Append:end);
                trackedFeatureNum(trackAdded,time2Append:end) = 0;

                %modify the matrix indicating linked feature information
                trackedFeatureInfo(track2Append,6*(time2Append-1)+1:end) = ...
                    trackedFeatureInfo(trackAdded,6*(time2Append-1)+1:end);
                trackedFeatureInfo(trackAdded,6*(time2Append-1)+1:end) = NaN;

            elseif mergeSplit && link21(i) > n && link21(i) <= n + numSplit %if this start is a split

                %get the track it split from
                trackSplitFrom = indxSplit(link21(i)-n);

                %find the tracks that have previously merged with this
                %splitting track
                mergeTrackToLookAt = find(indxMerge==trackSplitFrom);
                trackPrevMerge = link21(mergeTrackToLookAt+m);
                trackPrevMerge = trackPrevMerge(find(trackPrevMerge<=n));
                numPrevMerge = length(trackPrevMerge);

                %only consider splits that have possible previous merges
                if numPrevMerge ~= 0

                    %get the time of splitting
                    timeSplit = trackStartTime(i);

                    %find the times of merging
                    timeMerge = trackEndTime(trackPrevMerge) + 1;

                    %collect splitting and merging tracks into one
                    %matrix
                    trackedFeatMS = trackedFeatInfoInit([indxStart(i);...
                        indxEnd(trackPrevMerge)],:);

                    %find the costs for linking the merging tracks
                    %to the splitting track
                    costMatParams = costMatrices(4).costMatParam;
                    costMatParams.trackStats = trackStats;
                    eval(['[costVec,errFlag] = ' costMatrices(4).costMatFun ...
                        '(trackedFeatMS,timeSplit,timeMerge,'...
                        'costMatParams,gapCloseParam);'])

                    %choose the most likely merging track as that
                    %with the minimum cost
                    minCost = min(costVec);
                    mostProbTrack = find(costVec==minCost);
                    timeMerge = timeMerge(mostProbTrack);

                    if ~isinf(minCost)

                        %get the row number where this track is stored
                        %in trackedFeatInfoInit;
                        mostProbTrack = indxEnd(trackPrevMerge(mostProbTrack));

                        %modify the matrix indicating linked feature number
                        trackedFeatureNum(mostProbTrack,timeMerge:end) = ...
                            [trackedFeatureNum(trackSplitFrom,timeMerge:timeSplit-1) ...
                            trackedFeatureNum(indxStart(i),timeSplit:end)];
                        trackedFeatureNum(indxStart(i),timeSplit:end) = 0;

                        %modify the matrix indicating linked feature information
                        trackedFeatureInfo(mostProbTrack,6*(timeMerge-1)+1:end) = ...
                            [trackedFeatureInfo(trackSplitFrom,6*(timeMerge-1)+1:6*(timeSplit-1)) ...
                            trackedFeatureInfo(indxStart(i),6*(timeSplit-1)+1:end)];
                        trackedFeatureInfo(indxStart(i),6*(timeSplit-1)+1:end) = NaN;

                    end %(if ~isinf(minCost))

                end %(if numPrevMerge ~= 0)

            end %(if link21 <= n ... elseif ... mergeSplit && link21 ...)
            
        end %(for i=m:-1:1)

        %remove rows that do not contain tracks
        indx = find(max(trackedFeatureNum,[],2)>0);
        trackedFeatureNum = trackedFeatureNum(indx,:);
        trackedFeatureInfo = trackedFeatureInfo(indx,:);

    else %if there aren't any gaps to close

        trackedFeatureNum = trackedFeatNumInit;
        trackedFeatureInfo = trackedFeatInfoInit;

    end %(if n~=0 && m~=0 ... else ...)

    %get track statistics after this iteration of tracking
    eval(['[trackStatsT,statsRelChange,errFlag] = ' trackStatFun ...
        '(trackedFeatureInfo,lenFrac,timeWindow,trackStats);'])

    if errFlag %exit at this point if statistical analysis failed
        disp('--trackWithGapClosing: Getting track statistics failed. Stopping prematurely!');
        iterate = 0;
    elseif statsRelChange < tolerance %stop iterating if change is smaller than tolerance
        iterate = 0;
    end

end %(while iterate)


%%%%% ~~ the end ~~ %%%%%
