function [tracksFeatIndxCG,tracksCoordAmpCG,tracksFeatIndxLink,...
    tracksCoordAmpLink,kalmanInfoLink,errFlag] = trackCloseGapsKalman(...
    movieInfo,costMatParam,gapCloseParam,kalmanInitParam,useLocalDensity,...
    saveResults)
%TRACKCLOSEGAPSKALMAN links features between frames and closes gaps, with merging and splitting, using the Kalman filter
%
%SYNOPSIS [tracksFeatIndxCG,tracksCoordAmpCG,tracksFeatIndxLink,...
%    tracksCoordAmpLink,kalmanInfoLink,errFlag] = trackCloseGapsKalman(...
%    movieInfo,costMatParam,gapCloseParam,kalmanInitParam,useLocalDensity,...
%    saveResults)
%
%INPUT  movieInfo    : Array of size equal to the number of time points
%                      in a movie, containing the fields:
%             .xCoord      : Image coordinate system x-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%             .yCoord      : Image coordinate system y-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D. Default: zeros.
%             .zCoord      : Image coordinate system z-coordinates of detected
%                            features (in pixels). 1st column for
%                            value and 2nd column for standard deviation.
%                            Optional. Skipped if problem is 1D or 2D. Default: zeros.
%             .amp         : Amplitudes of PSFs fitting detected features.
%                            1st column for values and 2nd column
%                            for standard deviations.
%       costMatParam : Parameters needed for cost matrices. See cost
%                      matrices for details of fields.
%       gapCloseParam: Structure containing variables needed for gap closing.
%                      Contains the fields:
%             .timeWindow   : Largest time gap between the end of a track and the
%                             beginning of another that could be connected to it.
%             .tolerance    : Relative change in number of tracks in two
%                             consecutive gap closing steps below which
%                             iteration stops.
%             .mergeSplit   : Logical variable with value 1 if the merging
%                             and splitting of trajectories are to be consided;
%                             and 0 if merging and splitting are not allowed.
%       kalmanInitParam: Structure with fields containing variables
%                        used in Kalman filter initialization. See
%                        particular initialization function for fields.
%                        Optional. Enter [] or nothing if not to be used.
%       useLocalDensity: Structure with fields:
%           .link         : 1 if local density is used for determining
%                           search radius when linking between frames, 0
%                           otherwise.
%           .cg           : 1 if local density is used for determining
%                           search radius when gap closing, 0 otherwise.
%                        Both optional. Default: 0.
%       saveResults  : Structure with fields:
%           .dir          : Directory where results should be saved.
%                           Optional. Default: current directory.
%           .filename     : Name of file where results should be saved.
%                           Optional. Default: trackedFeatures.
%                       Whole structure optional.
%
%       All optional variables can be entered as [] to use default values.
%
%
%OUTPUT tracksFeatIndxCG  : Connectivity matrix of features between time points after gap closing.
%                           Rows indicate continuous tracks, while columns
%                           indicate time points. A track that ends before the
%                           last time point is followed by zeros, and a track
%                           that starts at a time after the first time point
%                           is preceded by zeros.
%       tracksCoordAmpCG  : The positions and amplitudes of the tracked
%                           features after gap closing. Number of rows = number of tracks,
%                           while number of columns = 8*number of time points.
%                           Each row consists of
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points
%                           where the track does not exist.
%       tracksFeatIndxLink: Connectivity matrix of features between time points after initial linking.
%                           Rows indicate continuous tracks, while columns
%                           indicate time points. A track that ends before the
%                           last time point is followed by zeros, and a track
%                           that starts at a time after the first time point
%                           is preceded by zeros.
%       tracksCoordAmpLink: The positions and amplitudes of the tracked
%                           features after initial linking. Number of rows = number of tracks,
%                           while number of columns = 8*number of time points.
%                           Each row consists of
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points
%                           where the track does not exist.
%       kalmanInfoLink    : Structure array with number of entries equal to
%                           number of frames in movie. Contains the fields:
%             .stateVec        : Kalman filter state vector for each
%                                feature in frame.
%             .stateCov        : Kalman filter state covariance matrix
%                                for each feature in frame.
%             .noiseVar        : Variance of state noise for each
%                                feature in frame.
%             .stateNoise      : Estimated state noise for each feature in
%                                frame.
%             .scheme          : 1st column: propagation scheme connecting
%                                feature to previous feature. 2nd column:
%                                propagation scheme connecting feature to
%                                next feature.
%       errFlag           : 0 if function executes normally, 1 otherwise.
%
%REMARKS The algorithm can handle cases where some frames do not have any
%        features at all. However, the very first frame must have some
%        features in it.
%
%Khuloud Jaqaman, April 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tracksFeatIndxCG   = [];
tracksCoordAmpCG   = [];
tracksFeatIndxLink = [];
tracksCoordAmpLink = [];
kalmanInfoLink     = [];
errFlag            = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 3
    disp('--trackCloseGapsKalman: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%get number of frames in movie
numFrames = length(movieInfo);

%check whether additional parameters for Kalman filter initialization are
%supplied
if nargin < 4 || isempty(kalmanInitParam)
    kalmanInitParam = [];
end

%determine whether local density is used
if nargin < 5 || isempty(useLocalDensity)
    useLocalDensity.link = 0;
    useLocalDensity.cg = 0;
else
    if ~isfield(useLocalDensity,'link')
        useLocalDensity.link = 0;
    end
    if ~isfield(useLocalDensity,'cg')
        useLocalDensity.cg = 0;
    end
end

%determine where to save results
if nargin < 6 || isempty(saveResults) %if nothing was input
    saveResDir = pwd;
    saveResFile = 'trackedFeatures';
else
    if ~isfield(saveResults,'dir') || isempty(saveResults.dir)
        saveResDir = pwd;
    else
        saveResDir = saveResults.dir;
    end
    if ~isfield(saveResults,'filename') || isempty(saveResults.filename)
        saveResFile = 'trackedFeatures';
    else
        saveResFile = saveResults.filename;
    end
end

%get gap closing parameters from input
tolerance = gapCloseParam.tolerance;

%calculate nearest neighbor distances and store them in movieInfo
for iFrame = 1 : numFrames

    %collect feature coordinates into one matrix
    coordinates = [movieInfo(iFrame).xCoord(:,1) movieInfo(iFrame).yCoord(:,1)];

    %compute distance matrix
    nnDist = createDistanceMatrix(coordinates,coordinates);

    %sort distance matrix and find nearest neighbor distance
    nnDist = sort(nnDist,2);
    nnDist = nnDist(:,2);

    %store nearest neighbor distance
    movieInfo(iFrame).nnDist = nnDist;

end %(for iFrame = 1 : numFrames)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Link between frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get initial tracks by linking features between consecutive frames
[tracksFeatIndxLink,tracksCoordAmpLink,kalmanInfoLink,nnDistLinkedFeat,...
    errFlag] = linkFeaturesKalman(movieInfo,costMatParam,[],...
    kalmanInitParam,useLocalDensity.link,useLocalDensity.nnWindowL);

%redo the linking by going backwards in the movie and using the
%Kalman filter information from the first linking attempt
%this will improve the linking and the state estimation
[tracksFeatIndxLink,tracksCoordAmpLink,kalmanInfoLink,nnDistLinkedFeat,...
    errFlag] = linkFeaturesKalman(movieInfo(end:-1:1),costMatParam,...
    kalmanInfoLink(end:-1:1),kalmanInitParam,useLocalDensity.link,...
    useLocalDensity.nnWindowL);

%go forward one more time to get the final estimate of the initial tracks
[tracksFeatIndxLink,tracksCoordAmpLink,kalmanInfoLink,nnDistLinkedFeat,...
    errFlag] = linkFeaturesKalman(movieInfo,costMatParam,...
    kalmanInfoLink(end:-1:1),kalmanInitParam,useLocalDensity.link,...
    useLocalDensity.nnWindowL);

%get number of tracks
numTracksOld = size(tracksFeatIndxLink,1);

%construct initial track connectivity matrix (which actually has no
%connectivity after the initial linking)
trackConnectLink = zeros(numTracksOld,numFrames);
trackConnectLink(:,1) = (1:numTracksOld)';

%find track start and end times
trackSEL = getTrackSEL(tracksCoordAmpLink);
trackStartTime = trackSEL(:,1);
trackEndTime = trackSEL(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Close gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the indices and total number of tracks that will be considered for gap
%closing (those which start after the first frame)
indxStart = find(trackStartTime > 1);
numStart = length(indxStart);

%get the indices and total number of tracks that will be considered for gap
%closing (those which end before the last frame)
indxEnd = find(trackEndTime < numFrames);
numEnd = length(indxEnd);

%initialize track connectivity matrix after gap closing = track
%connectivity matrix from initial linking
trackConnectCGOld = trackConnectLink;

%if there are gaps to close ...
if numEnd ~= 0 && numStart ~= 0

    %assign dummy value to relChangeNumTracks to start first iteration
    relChangeNumTracks = 1.5;

    while relChangeNumTracks > tolerance

        %calculate the cost matrix, which already includes the
        %costs of birth and death
        [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
            errFlag] = costMatLinearMotionCloseGaps(tracksCoordAmpLink,...
            tracksFeatIndxLink,indxStart,indxEnd,trackStartTime,...
            trackEndTime,costMatParam,gapCloseParam,kalmanInfoLink,...
            trackConnectCGOld,useLocalDensity.cg,nnDistLinkedFeat,...
            useLocalDensity.nnWindowCG);

        %link tracks based on this cost matrix, allowing for birth and death
        [link12,link21] = lap(costMat,nonlinkMarker);

        %initialize new track connectivity matrix
        trackConnectCGNew = trackConnectLink;

        for iStart = numStart : -1 : 1 %go over all track starts

            if link21(iStart) <= numEnd %if this start is linked to an end

                %get the track to be appended at its end, the time point at which
                %it will be appended, and the track that will be added to it
                track2Append = indxEnd(link21(iStart));
                trackAdded = indxStart(iStart);

                %update track connectivity matrix
                trackConnectCGNew(track2Append,2:end) = trackConnectCGNew(trackAdded,1:end-1);
                trackConnectCGNew(trackAdded,:) = 0;

            end %(if link21(iStart) <= numEnd)

        end %(for iStart = numStart : -1 : 1)

        %remove rows that do not contain tracks
        maxValue = max(trackConnectCGNew,[],2);
        trackConnectCGNew = trackConnectCGNew(maxValue > 0,:);

        %determine number of tracks after gap closing
        numTracksNew = size(trackConnectCGNew,1);

        %calculate relative change in number of tracks from previous
        %iteration
        relChangeNumTracks = abs(numTracksNew - numTracksOld) / numTracksOld;
        disp([num2str(numTracksNew) '   ' num2str(relChangeNumTracks)]);
        
        %change New to Old
        trackConnectCGOld = trackConnectCGNew;
        numTracksOld = numTracksNew;

    end %(while relChangeNumTracks > tolerance)

    %reserve memory for matrix storing feature indices
    %and for matrix storing coordinates and amplitudes
    tracksFeatIndxCG = zeros(numTracksNew,numFrames);
    tracksCoordAmpCG = NaN*ones(numTracksNew,8*numFrames);

    %go over all tracks
    for iTrack = 1 : numTracksNew

        %get indices of segments from linking making up current track
        segmentIndx = trackConnectCGNew(iTrack,:);
        segmentIndx = segmentIndx(segmentIndx ~= 0);
        
        %go over all segments
        for iSegment = segmentIndx

            %get segment start and end time
            segmentStartTime = trackStartTime(iSegment);
            segmentEndTime = trackEndTime(iSegment);
            
            %store feature indices
            tracksFeatIndxCG(iTrack,segmentStartTime:segmentEndTime) = ...
                tracksFeatIndxLink(iSegment,segmentStartTime:segmentEndTime);

            %modify start and end time for storing coordinates
            %and amplitudes
            segmentStartTime = 8 * (segmentStartTime - 1) + 1;
            segmentEndTime = 8 * segmentEndTime;

            %store coordinates and amplitudes
            tracksCoordAmpCG(iTrack,segmentStartTime:segmentEndTime) = ...
                tracksCoordAmpLink(iSegment,segmentStartTime:segmentEndTime);

        end %(for iSegment = segmentIndx)
        
    end %(for iTrack = 1 : numTracksNew)

end %(if n~=0 && m~=0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([saveResDir filesep saveResFile],'costMatParam','gapCloseParam',...
    'kalmanInitParam','tracksFeatIndxCG','tracksCoordAmpCG',...
    'tracksFeatIndxLink','tracksCoordAmpLink','kalmanInfoLink');


%%%%% ~~ the end ~~ %%%%%


% %
% %
% %
% %             elseif mergeSplit && link21(i) > n && link21(i) <= n + numSplit %if this start is a split
% %
% %                 %get the track it split from
% %                 trackSplitFrom = indxSplit(link21(i)-n);
% %
% %                 %find the tracks that have previously merged with this
% %                 %splitting track
% %                 mergeTrackToLookAt = find(indxMerge==trackSplitFrom);
% %                 trackPrevMerge = link21(mergeTrackToLookAt+m);
% %                 trackPrevMerge = trackPrevMerge(trackPrevMerge<=n);
% %                 numPrevMerge = length(trackPrevMerge);
% %
% %                 %only consider splits that have possible previous merges
% %                 if numPrevMerge ~= 0
% %
% %                     %get the time of splitting
% %                     timeSplit = trackStartTime(i);
% %
% %                     %find the times of merging
% %                     timeMerge = trackEndTime(trackPrevMerge) + 1;
% %
% %                     %collect splitting and merging tracks into one
% %                     %matrix
% %                     trackedFeatMS = trackedFeatureInfo([indxStart(i);...
% %                         indxEnd(trackPrevMerge)],:);
% %
% %                     %find the costs for linking the merging tracks
% %                     %to the splitting track
% %                     costMatParams = costMatrices(4).costMatParam;
% %                     costMatParams.trackStats = trackStats;
% %                     eval(['[costVec,errFlag] = ' costMatrices(4).costMatFun ...
% %                         '(trackedFeatMS,timeSplit,timeMerge,'...
% %                         'costMatParams,gapCloseParam);'])
% % 
% %                     %choose the most likely merging track as that
% %                     %with the minimum cost
% %                     minCost = min(costVec);
% %                     mostProbTrack = find(costVec==minCost);
% %                     timeMerge = timeMerge(mostProbTrack);
% % 
% %                     if ~isinf(minCost)
% % 
% %                         %get the row number where this track is stored
% %                         %in trackedFeatureInfo;
% %                         mostProbTrack = indxEnd(trackPrevMerge(mostProbTrack));
% % 
% %                         %modify the matrix indicating linked feature number
% %                         trackedFeatureIndx(mostProbTrack,timeMerge:end) = ...
% %                             [trackedFeatureIndx(trackSplitFrom,timeMerge:timeSplit-1) ...
% %                             trackedFeatureIndx(indxStart(i),timeSplit:end)];
% %                         trackedFeatureIndx(indxStart(i),timeSplit:end) = 0;
% % 
% %                         %modify the matrix indicating linked feature information
% %                         trackedFeatureInfo(mostProbTrack,8*(timeMerge-1)+1:end) = ...
% %                             [trackedFeatureInfo(trackSplitFrom,8*(timeMerge-1)+1:8*(timeSplit-1)) ...
% %                             trackedFeatureInfo(indxStart(i),8*(timeSplit-1)+1:end)];
% %                         trackedFeatureInfo(indxStart(i),8*(timeSplit-1)+1:end) = NaN;
% % 
% %                     end %(if ~isinf(minCost))
% % 
% %                 end %(if numPrevMerge ~= 0)
% % 


% %close gaps sequentially (in segments) to avoid memory issues if problem is too big
% %determine the lower and upper bounds of first segment for gap closing
% %go back in time by timeWindow when looking for track ends to be
% %connected to track starts in this segment
% segmentUB = numFrames; %upper bound (for both ends and starts)
% segmentLBS = max(segmentUB+1-segmentLength,2); %lower bound for starts
% segmentLBE = max(segmentUB+1-segmentLength-timeWindow,1); %lower bounds for ends
% 
% %if there is no merging and splitting, segmentUB cannot be smaller than
% %3, because if segmentUB < 3, there won't be any gaps to close
% %if there is merging and splitting, segmentUB cannot be smaller than 2,
% %because if segmentUB < 2, there won't be any gaps to close or merges
% %and splits to consider
% while segmentUB > 2 || (mergeSplit && segmentUB > 1)
% 
%     %get the index of tracks that start between segment start time and
%     %segment end time
%     trackStartTime = trackSEL(:,1);
%     indxStart = find(trackStartTime >= segmentLBS & trackStartTime <= segmentUB);
%     numStart = length(indxStart);
% 
%     %get the index of tracks that end between segment start time -
%     %timeWindow and segment end time
%     trackEndTime = trackSEL(:,2);
%     indxEnd = find(trackEndTime >= segmentLBE & trackEndTime <= min(segmentUB,numFrames-1));
%     numEnd = length(indxEnd);
% 
%     %go to next segment. Update lower and upper bounds
%     segmentUB = segmentLBS - 1; %upper bound (for both end and starts)
%     segmentLBS = max(segmentUB+1-segmentLength,1); %lower bound for starts
%     segmentLBE = max(segmentUB+1-segmentLength-timeWindow,1); %lower bounds for ends
% 
% end %(segmentUB > 2 || (mergeSplit && segmentUB > 1))

%get parameters from input
% timeWindow = gapCloseParam.timeWindow;
% mergeSplit = gapCloseParam.mergeSplit;
% if isfield(gapCloseParam,'segmentLength')
%     segmentLength = gapCloseParam.segmentLength;
% else
%     segmentLength = numFrames;
% end

