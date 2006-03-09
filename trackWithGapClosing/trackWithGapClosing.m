function [trackedFeatureNum,trackedFeatureInfo,errFlag] = trackWithGapClosing(...
    movieInfo,linkingLimits,gapCloseParam,iterParam)
%TRACKWITHGAPCLOSING links features between frames in a movie with gap closing, with the possibility of track merging and splitting
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
%       linkingLimits: Structure with variables defining limits beyond which
%                      linking is not possible. Contains the fields:
%             .searchRadius: Maximum distance between two features in two
%                            consecutive time points that allows linking 
%                            them (in pixels) in initial linking.
%             .maxAmpRatio : Maximum ratio between the amplitudes of two
%                            features in two censecutive time points that 
%                            allows linking them in initial linking.
%             .cutoffCProb : Cumulative probability of squared displacement
%                            or amplitude change beyond which links are not
%                            allowed in linking based on statistical
%                            analysis of tracks.
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
if nargin ~= 4
    disp('--trackWithGapClosing: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

searchRadius = linkingLimits.searchRadius;
maxAmpRatio  = linkingLimits.maxAmpRatio;
cutoffCProb  = linkingLimits.cutoffCProb;
timeWindow   = gapCloseParam.timeWindow;
mergeSplit   = gapCloseParam.mergeSplit;
tolerance    = iterParam.tolerance;
lenFrac      = iterParam.lenFrac;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of time points in movie
numTimePoints = length(movieInfo);

%get initial tracks by linking features between consecutive frames using
%the simple initial criteria
costMatParams = struct('searchRadius',searchRadius,'maxAmpRatio',maxAmpRatio);
[trackedFeatNumInit,trackedFeatInfoInit,errFlag] = linkFeaturesTp2Tp(...
    movieInfo,'costMatSimple',costMatParams);

%get track statistics
[dispSqLambdaT,ampDiffStdT,errFlag] = getTrackStats(...
    trackedFeatInfoInit,lenFrac,timeWindow);

%exit at this point if statistical analysis could not be performed
if errFlag
    disp('--trackWithGapClosing: getTrackStats failed. Stopping after initial linking!');
    iterate = 0;
else %otherwise iterate
    iterate = 1;
end

while iterate

    %get latest statistics
    dispSqLambda = dispSqLambdaT;
    ampDiffStd = ampDiffStdT;

    %get the maximum squared displacement that allows linking 2 features
    maxDispSq = expinv(cutoffCProb,dispSqLambda);

    %find the maximum squared amplitude change that allows linking 2 features
    maxAmpDiffSq = (norminv(cutoffCProb,0,ampDiffStd)).^2;

    %calculate the additive constant for each cost
    addConst = log(ampDiffStd) - log(dispSqLambda);    
    
    %link features between consecutive frames using the information obtained
    %via the statistical analysis of tracks
    costMatParams = struct('dispSqLambda',dispSqLambda(1),'ampDiffStd',...
        ampDiffStd(1),'maxDispSq',maxDispSq(1),'maxAmpDiffSq',maxAmpDiffSq(1));
    [trackedFeatNumInit,trackedFeatInfoInit,errFlag] = linkFeaturesTp2Tp(...
        movieInfo,'costMatLogL',costMatParams);

    %close gaps using the information obtained via the statistical analysis
    %of tracks ...

    %get number of tracks formed by initial linking
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
                    dispSq = (xCoordEnd(i)-xCoordStart(j))^2 + ...
                        (yCoordEnd(i)-yCoordStart(j))^2;

                    %calculate the square difference between the amplitude at
                    %the beginning of track j and the end of track i
                    ampDiffSq = (ampEnd(i) - ampStart(j))^2;

                    %if this is a possible link ...
                    if dispSq < maxDispSq(timeGap) && ...
                            ampDiffSq < maxAmpDiffSq(timeGap)

                        %assign the cost for this pair of tracks
                        indx1 = [indx1; i]; %row number
                        indx2 = [indx2; j]; %column number
                        cost = [cost; dispSqLambda(timeGap)*dispSq + ...
                            ampDiffSq/2/ampDiffStd(timeGap)^2 + ...
                            addConst(timeGap)]; %cost

                    end %(if dispSq < maxDispSq(timeGap) && ampDIffSq < maxAmpDiffSq(timeGap))

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
        [link12,link21] = lap(costMat,-1000,0,1);

        %Connect the tracks using the information from LAP
        trackedFeatureNum  = trackedFeatNumInit;
        trackedFeatureInfo = trackedFeatInfoInit;

        %if merging and splitting are considered ...
        if mergeSplit

            for i=m:-1:1 %go over all track starts
                if link21(i) > n && link21(i) <= n + numSplit %if this start is a split

                    %get the track it split from
                    trackSplitFrom = indxSplit(link21(i)-n);

                    %find the tracks that have previously merged with this
                    %splitting track
                    mergeTrackToLookAt = find(indxMerge==trackSplitFrom);
                    trackPrevMerge = link21(mergeTrackToLookAt+m);
                    trackPrevMerge = trackPrevMerge(find(trackPrevMerge<=n));

                    %only consider splits that have possible previous merges
                    if ~isempty(trackPrevMerge)

                        %get the time of splitting
                        timeSplit = trackStartTime(i);

                        %find the times of merging
                        timeMerge = trackEndTime(trackPrevMerge)+1;

                        %calculate the ratio of the amplitudes before merging
                        %and the amplitude after splitting (max/min)
                        ampRatio = ampEnd(trackPrevMerge)/ampStart(i);
                        ampRatio = max([ampRatio 1./ampRatio])';

                        %retain only those merges which happened before the
                        %splitting and whose amplitudes at the end are close
                        %to the amplitude of the splitting track
                        tmp = find(timeMerge<timeSplit & ampRatio<maxAmpRatioM);
                        timeMerge = timeMerge(tmp);
                        trackPrevMerge = trackPrevMerge(tmp);

                        %if there really are merges before this split with
                        %close enough amplitude...
                        if ~isempty(trackPrevMerge)

                            %choose the track whose merging time is closest to
                            %this track's splitting time
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

        %close gaps due to out-of-focus disappearance
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

    %get track statistics after this iteration of tracking
    [dispSqLambdaT,ampDiffStdT,errFlag] = getTrackStats(...
        trackedFeatureInfo,lenFrac,timeWindow);

    if errFlag %exit at this point if statistical analysis failed

        disp('--trackWithGapClosing: getTrackStats failed. Stopping prematurely!');
        iterate = 0;

    else %otherwise check change in parameters

        relParamChange = [abs((dispSqLambdaT-dispSqLambda)./dispSqLambda); ...
            abs((ampDiffStdT-ampDiffStd)./ampDiffStd)];

        %stop iterating if change is smaller than tolerance
        if min(relParamChange < tolerance)
            iterate = 0;
        end

    end %(if errFlag ... else ...)


end %(while iterate)


%%%%% ~~ the end ~~ %%%%%
