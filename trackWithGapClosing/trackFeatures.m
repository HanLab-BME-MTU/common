function [trackedFeatureNum,trackedFeatureInfo,errFlag] = trackFeatures(...
    movieInfo,tp2tpLinkCriteria,gapClosingCriteria)
%TRACKFEATURES links features between frames in a movie with gap closing, with the possibility of track merging and splitting
%
%SYNOPSIS [trackedFeatureNum,trackedFeatureInfo,errFlag] = trackFeatures(...
%    movieInfo,tp2tpLinkCriteria,gapClosingCriteria)
%
%INPUT  movieInfo        : Array of size equal to the number of time points
%                          in a movie, containing the fields:
%             .xCoord       : Image coordinate system x-coordinate of detected
%                             features (in pixels).
%             .yCoord       : Image coorsinate system y-coordinate of detected
%                             features (in pixels).
%             .amp          : Amplitudes of PSFs fitting detected features.
%       tp2tpLinkCriteria: Structure with variables related to initial
%                          linking of features between frames. Contains the
%                          fields:
%             .searchRadius1: Maximum distance between two features in two
%                             consecutive time points that allows linking 
%                             them (in pixels).
%             .maxAmpRatio1 : Maximum ratio between the amplitudes of two
%                             features in two censecutive time points that 
%                             allows linking them.
%       gapClosingCriteria:Structure containing additional variables that
%                          are needed for gap closing. 
%                          Optional. Use [] for default.
%                          Contains the fields:
%             .timeWindow   : Largest time gap between the end of a track and the
%                             beginning of another that could be connected to it.
%                             Default: 10. Can skip if default.
%             .mergeSplit   : Logical variable with value 1 if the merging 
%                             and splitting of trajectories are to be consided;
%                             and 0 if merging and splitting are not allowed.
%                             Default: 0. Can skip if default.
%             .searchRadiusM: Maximum search radius, relevant when motion is confined
%                             within a volume (in pixels). 
%                             Default: infinity. Can skip if default.
%             .maxAmpRatioM : Maximum ratio between the amplitudes of two 
%                             features separated by "timeWindow" time
%                             points apart that allows connecting them. 
%                             Default: maxAmpRatio1. Can skip if default.
%             .minAmpRatioMS: Minimum acceptable ratio between the
%                             amplitude of the superposition of two 
%                             features and the sum of the amplitudes of 
%                             the two features before merging or after
%                             splitting.
%                             Default: 0.9. Can skip if default.
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
%Khuloud Jaqaman, February 2006

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
if nargin < 2
    disp('--trackFeatures: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check tp2tpLinkCriteria
if ~isfield(tp2tpLinkCriteria,'searchRadius1')
    disp('--trackFeatures: Variable tp2tpLinkCriteria should contain the field searchRadius1!');
    errFlag = 1;
else
    searchRadius1 = tp2tpLinkCriteria.searchRadius1;
    if searchRadius1 <= 0
        disp('--trackFeatures: tp2tpLinkCriteria.searchRadius1 should be positive!');
        errFlag = 1;
    end
end
if ~isfield(tp2tpLinkCriteria,'maxAmpRatio1')
    disp('--trackFeatures: Variable tp2tpLinkCriteria should contain the field maxAmpRatio1');
    errFlag = 1;
else
    maxAmpRatio1 = tp2tpLinkCriteria.maxAmpRatio1;
    if maxAmpRatio1 <= 1
        disp('--trackFeatures: tp2tpLinkCriteria.maxAmpRatio1 should be > 1!');
        errFlag = 1;
    end
end

%exit if there are problems with input
if errFlag
    disp('--trackFeatures: Please fix input data!')
    return
end

%assign defaults
timeWindow_def = 10;
mergeSplit_def = 0;
searchRadiusM_def = Inf;
maxAmpRatioM_def = maxAmpRatio1;
minAmpRatioMS_def = 0.9;

%check gapClosingCriteria
if nargin < 3 || isempty(gapClosingCriteria)
    
    timeWindow    = timeWindow_def;
    mergeSplit    = mergeSplit_def;
    searchRadiusM = searchRadiusM_def;
    maxAmpRatioM  = maxAmpRatioM_def;
    minAmpRatioMS = minAmpRatioMS_def;

else

    if ~isfield(gapClosingCriteria,'timeWindow') || ...
            isempty(gapClosingCriteria.timeWindow)
        timeWindow = timeWindow_def;
    else
        timeWindow = gapClosingCriteria.timeWindow;
        if timeWindow <= 1
            disp('--trackFeatures: Variable "timeWindow" should be larger than 1!');
            errFlag = 1;
        end
    end

    if ~isfield(gapClosingCriteria,'mergeSplit') || ...
            isempty(gapClosingCriteria.mergeSplit)
        mergeSplit = mergeSplit_def;
    else
        mergeSplit = gapClosingCriteria.mergeSplit;
        if mergeSplit~=0 && mergeSplit~=1
            disp('--trackFeatures: Variable "mergeSplit" should be either 0 or 1!');
            errFlag = 1;
        end
    end

    if ~isfield(gapClosingCriteria,'searchRadiusM') || ...
            isempty(gapClosingCriteria.searchRadiusM)
        searchRadiusM = searchRadiusM_def;
    else
        searchRadiusM = gapClosingCriteria.searchRadiusM;
        if searchRadiusM <= 0
            disp('--trackFeatures: gapClosingCriteria.searchRadiusM should be positive!');
            errFlag = 1;
        end
    end
    
    if ~isfield(gapClosingCriteria,'maxAmpRatioM') || ...
            isempty(gapClosingCriteria.maxAmpRatioM)
        maxAmpRatioM = maxAmpRatioM_def;
    else
        maxAmpRatioM = gapClosingCriteria.maxAmpRatioM;
        if maxAmpRatioM <= 1
            disp('--trackFeatures: gapClosingCriteria.maxAmpRatioM should be > 1!');
            errFlag = 1;
        end
    end
    
    if ~isfield(gapClosingCriteria,'minAmpRatioMS') || ...
            isempty(gapClosingCriteria.minAmpRatioMS)
        minAmpRatioMS = minAmpRatioMS_def;
    else
        minAmpRatioMS = gapClosingCriteria.minAmpRatioMS;
        if minAmpRatioMS > 1 || minAmpRatioMS <= 0
            disp('--trackFeatures: gapClosingCriteria.minAmpRatioMS should be > 0 & < 1!');
            errFlag = 1;
        end
    end
        
end

%exit if there are problems with input
if errFlag
    disp('--trackFeatures: Please fix input data!')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get initial tracks by linking features between consecutive frames
%note that these tracks are arranged in ascending order of the time they
%start and, within each group starting at the same time point, in
%descending order of their length
[trackedFeatNumInit,trackedFeatInfoInit,errFlag] = linkFeaturesTp2Tp( ...
    movieInfo,searchRadius1,maxAmpRatio1);

%get number of time points in movie
numTimePoints = length(movieInfo);

%get number of tracks
numTracks = size(trackedFeatNumInit,1);

%find the starting time points of all tracks
trackStartTime = zeros(numTracks,1);
for i=1:numTracks
    trackStartTime(i) = find((trackedFeatNumInit(i,:)~=0),1,'first');
end

%find the termination time points of all tracks
trackEndTime = zeros(numTracks,1);
for i=1:numTracks
    trackEndTime(i) = find((trackedFeatNumInit(i,:)~=0),1,'last');
end

%get the x,y-coordinates and amplitudes of features at the starts of
%their tracks
xCoordStart = zeros(numTracks,1);
yCoordStart = zeros(numTracks,1);
ampStart = zeros(numTracks,1);
for i=1:numTracks
    xCoordStart(i) = trackedFeatInfoInit(i,(trackStartTime(i)-1)*6+1);
    yCoordStart(i) = trackedFeatInfoInit(i,(trackStartTime(i)-1)*6+2);
    ampStart(i) = trackedFeatInfoInit(i,(trackStartTime(i)-1)*6+3);
end

%get the x,y-coordinates and amplitudes of features at the ends of
%their tracks
xCoordEnd = zeros(numTracks,1);
yCoordEnd = zeros(numTracks,1);
ampEnd = zeros(numTracks,1);
for i=1:numTracks
    xCoordEnd(i) = trackedFeatInfoInit(i,(trackEndTime(i)-1)*6+1);
    yCoordEnd(i) = trackedFeatInfoInit(i,(trackEndTime(i)-1)*6+2);
    ampEnd(i) = trackedFeatInfoInit(i,(trackEndTime(i)-1)*6+3);
end

%remove tracks that start at the first time point and get the total number of
%tracks whose starts are to be considered
indxStart = find(trackStartTime>1);
trackStartTime = trackStartTime(indxStart);
xCoordStart = xCoordStart(indxStart);
yCoordStart = yCoordStart(indxStart);
ampStart = ampStart(indxStart);
m = length(indxStart);

%remove tracks that end at the last time point and get the total number of
%tracks whose ends are to be considered
indxEnd = find(trackEndTime<numTimePoints);
trackEndTime = trackEndTime(indxEnd);
xCoordEnd = xCoordEnd(indxEnd);
yCoordEnd = yCoordEnd(indxEnd);
ampEnd = ampEnd(indxEnd);
n = length(indxEnd);

%if there are gaps to close ...
if n~=0 && m~=0

    %get the square of the search radius as a funtion of time between
    %frames
    searchRadius1Sq = searchRadius1^2;
    searchRadiusMSq = searchRadiusM^2;
    for i=timeWindow:-1:1
        searchRadiusSq(i) = min(i*searchRadius1Sq,searchRadiusMSq);
    end
    
    %get the maximum allowed amplitude ratio as a function of time
    %between frames
    for i=timeWindow:-1:1
        maxAmpRatio(i) = maxAmpRatio1 + (i-1)*(maxAmpRatioM-maxAmpRatio1)...
            /(timeWindow-1);
    end

    %generate cost matrix, in sparse matrix format

    indx1 = []; %row number in cost matrix
    indx2 = []; %column number in cost matrix
    cost  = []; %cost value

    %costs for closing gaps due to features going out-of focus

    for j=1:m %go over all starts
        for i=1:n %go over all ends

            %subtract starting time from ending time
            timeGap = trackStartTime(j) - trackEndTime(i);

            %if the time gap between the end of track i and the start of track
            %j is within the acceptable limits ...
            if timeGap > 1 && timeGap <= timeWindow

                %calculate the square distance between the end
                %point and starting point
                distance2 = ((xCoordEnd(i)-xCoordStart(j))^2 + ...
                    (yCoordEnd(i)-yCoordStart(j))^2)/timeGap;

                %if the distance is smaller than the search 
                %radius relevent to that time gap ...
                if distance2 < searchRadiusSq(1)

                    %calculate the amplitude ratio
                    ampRatio = max(ampStart(j),ampEnd(i))/...
                        min(ampStart(j),ampEnd(i));

                    %if the amplitude ratio is smaller than the maximum
                    %allowed ratio for that time gap ...
                    if ampRatio < maxAmpRatio(timeGap)

                        %assign the cost for this pair of tracks
                        indx1 = [indx1; i]; %row number
                        indx2 = [indx2; j]; %column number
                        cost = [cost; distance2*ampRatio]; %cost
                        
                    end %(if ampRatio < maxAmpRatio(timeGap))

                end %(if distance2 < searchRadiusSq(timeGap))

            end %(if timeGap > 1 || timeGap <= timeWindow)

        end %(for i=1:n)
    end %(for j=1:m)

    %if merging and splitting are to be considered ...
    if mergeSplit

        %costs for closing gaps due to merging

        numMerge  =  0; %index counting merging events
        indxMerge = []; %vector storing merging track number

        for j=1:n %go over all ends
            for i=1:numTracks %go over all tracks

                timeIndx  = trackEndTime(j)*6;
                xCoordMid = trackedFeatInfoInit(i,timeIndx+1);
                yCoordMid = trackedFeatInfoInit(i,timeIndx+2);
                distance2 = (xCoordEnd(j)-xCoordMid)^2 ...
                    + (yCoordEnd(j)-yCoordMid)^2;

                %if the distance is smaller than the search radius ...
                if distance2 < searchRadiusSq(1)
                    
                    %calculate the ratio between the amplitude of the
                    %merged feature and the sum of the amplitudes of the
                    %merging ones
                    ampMidT0 = trackedFeatInfoInit(i,timeIndx-3); %amplitude before merging
                    ampMidT1 = trackedFeatInfoInit(i,timeIndx+3); %amplitude after merging
                    ampRatio = ampMidT1/(ampMidT0+ampEnd(j));
                    
                    %if this amplitude ratio is within reasonable limits ...
                    if ampRatio > minAmpRatioMS && ampRatio < maxAmpRatio1
                        
                        %increase the "merge index" by one
                        numMerge = numMerge + 1;

                        %save the merging track's number
                        indxMerge = [indxMerge; i];

                        %calculate the cost of merging
                        indx1 = [indx1; j]; %row number
                        indx2 = [indx2; numMerge+m]; %column number
                        cost = [cost; distance2*max(ampRatio,1/ampRatio)]; %cost
                        
                    end %(if ampRatio > 0.9 && ampRatio < maxAmpRatio1)

                end %(if distance2 < searchRadiusSq(1))

            end %(for i=1:numTracks)
        end %(for j=1:n)

        %costs for closing gaps due to splitting
        
        numSplit  =  0; %index counting splitting events
        indxSplit = []; %vector storing splitting track number

        for j=1:m %go over all starts
            for i=1:numTracks %go over all tracks

                timeIndx = (trackStartTime(j)-2)*6;
                xCoordMid = trackedFeatInfoInit(i,timeIndx+1);
                yCoordMid = trackedFeatInfoInit(i,timeIndx+2);
                distance2 = (xCoordStart(j)-xCoordMid)^2 ...
                    + (yCoordStart(j)-yCoordMid)^2;

                %if the distance is smaller than the search radius ...
                if distance2 < searchRadiusSq(1)

                    %calculate the ratio between the amplitude of the
                    %feature before splitting and the sum of the amplitudes
                    %of the features after splitting
                    ampMidT0 = trackedFeatInfoInit(i,timeIndx+3); %amplitude before splitting
                    ampMidT1 = trackedFeatInfoInit(i,timeIndx+9); %amplitude after splitting
                    ampRatio = ampMidT0/(ampMidT1+ampStart(j));

                    %if this amplitude ratio is within reasonable limits ...
                    if ampRatio > minAmpRatioMS && ampRatio < maxAmpRatio1

                        %increase the "split index" by one
                        numSplit = numSplit + 1;

                        %save the splitting track's number
                        indxSplit = [indxSplit; i];

                        %calculate the cost of splitting
                        indx1 = [indx1; numSplit+n]; %row number
                        indx2 = [indx2; j]; %column number
                        cost = [cost; distance2*max(ampRatio,1/ampRatio)]; %cost

                    end %(if ampRatio > 0.9 && ampRatio < maxAmpRatio1)

                end %(if distance2 < searchRadiusSq(1))

            end %(for i=1:numTracks)
        end %(for j=1:m)

    end %(if mergeSplit)

    %create cost matrix
    costMat = sparse(indx1,indx2,cost);

    %link tracks based on this cost matrix, allowing for birth and death
    [link12,link21] = lap(costMat,-1,0,1);

    %Connect the tracks using the information from LAP
    trackedFeatureNum = trackedFeatNumInit;
    trackedFeatureInfo = trackedFeatInfoInit;

    %if merging and splitting are to be considered ...
    if mergeSplit

        %connect splits and merges first
        for i=m:-1:1 %go over all track starts
            if link21(i) > n && link21(i) <= n + numSplit %if this start is a split

                %get the track it split from
                trackSplitFrom = indxSplit(link21(i)-n);

                %find the tracks that have previously merged with the
                %splitting track
                mergeTrackToLookAt = find(indxMerge==trackSplitFrom);
                trackPrevMerge = link21(mergeTrackToLookAt+m);
                trackPrevMerge = trackPrevMerge(find(trackPrevMerge<=n));

                %only consider splits that have a previous merge
                if ~isempty(trackPrevMerge)

                    %find the most likely track that had merged. Choose
                    %that whose merging time is closest to this track's
                    %splitting time
                    timeSplit = trackStartTime(i);
                    timeMerge = trackEndTime(trackPrevMerge)+1;
                    tmp = find(timeMerge<timeSplit);
                    timeMerge = timeMerge(tmp);
                    trackPrevMerge = trackPrevMerge(tmp);

                    %again, only consider splits that have a previous merge
                    if ~isempty(trackPrevMerge)

                        timeMerge = timeMerge(end);
                        trackPrevMerge = trackPrevMerge(end);

                        %modify the matrix indicating linked feature number
                        trackedFeatureNum(indxEnd(trackPrevMerge),timeMerge:end) = ...
                            [trackedFeatureNum(trackSplitFrom,timeMerge:timeSplit-1) ...
                            trackedFeatureNum(indxStart(i),timeSplit:end)];
                        trackedFeatureNum(indxStart(i),timeSplit:end) = 0;

                        %modify the matrix indicating linked feature information
                        trackedFeatureInfo(indxEnd(trackPrevMerge),6*(timeMerge-1)+1:end) = ...
                            [trackedFeatureInfo(trackSplitFrom,6*(timeMerge-1)+1:6*(timeSplit-1)) ...
                            trackedFeatureInfo(indxStart(i),6*(timeSplit-1)+1:end)];
                        trackedFeatureInfo(indxStart(i),6*(timeSplit-1)+1:end) = NaN;

                    end

                end

            end
        end

    end %(if mergeSplit)
    
    %then close gaps due to out-of-focus disappearance
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

        end
    end

    %remove rows that do not contain tracks
    indx = find(max(trackedFeatureNum,[],2)>0);
    trackedFeatureNum = trackedFeatureNum(indx,:);
    trackedFeatureInfo = trackedFeatureInfo(indx,:);

else %if there aren't any gaps to close

    trackedFeatureNum = trackedFeatNumInit;
    trackedFeatureInfo = trackedFeatInfoInit;

end %(if n~=0 && m~=0 ... else ...)


%%%%% ~~ the end ~~ %%%%%
