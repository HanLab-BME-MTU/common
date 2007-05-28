function [tracksFinal,tracksFeatIndxLink,...
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
%        When gapCloseParam.mergeSplit = 1, gapCloseParam.tolerance MUST BE
%        1 to prevent iteration. The part which assigns track segment types
%        based on overall track type cannot handle tracks which contain
%        merging and splitting. Therefore, it will only work well with the
%        initial connectivity matrix from linking.
%
%Khuloud Jaqaman, April 2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tracksFinal        = [];
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
% tolerance = gapCloseParam.tolerance;
mergeSplit = gapCloseParam.mergeSplit;

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
[tracksFeatIndxLink,tracksCoordAmpLink,kalmanInfoLink] = ...
    linkFeaturesKalman(movieInfo,costMatParam,[],...
    kalmanInitParam,useLocalDensity.link,useLocalDensity.nnWindowL);
clear tracksFeatIndxLink tracksCoordAmpLink errFlag

%redo the linking by going backwards in the movie and using the
%Kalman filter information from the first linking attempt
%this will improve the linking and the state estimation
[tracksFeatIndxLink,tracksCoordAmpLink,kalmanInfoLink] = ...
    linkFeaturesKalman(movieInfo(end:-1:1),costMatParam,...
    kalmanInfoLink(end:-1:1),kalmanInitParam,useLocalDensity.link,...
    useLocalDensity.nnWindowL);
clear tracksFeatIndxLink tracksCoordAmpLink errFlag

%go forward one more time to get the final estimate of the initial tracks
[tracksFeatIndxLink,tracksCoordAmpLink,kalmanInfoLink,nnDistLinkedFeat,...
    errFlag] = linkFeaturesKalman(movieInfo,costMatParam,...
    kalmanInfoLink(end:-1:1),kalmanInitParam,useLocalDensity.link,...
    useLocalDensity.nnWindowL);

%get number of tracks
numTracksLink = size(tracksFeatIndxLink,1);

%get track start and end times
trackSEL = getTrackSEL(tracksCoordAmpLink);
trackStartTime = trackSEL(:,1);
trackEndTime   = trackSEL(:,2);
clear trackSEL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Close gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if there are gaps to close (i.e. if there are tracks that start after the
%first frame and tracks that end before the last frame) ...
if any(trackStartTime > 1) && any(trackEndTime < numFrames)

    %calculate the cost matrix, which already includes the
    %costs of birth and death
    [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
        errFlag] = costMatLinearMotionCloseGaps(tracksCoordAmpLink,...
        tracksFeatIndxLink,trackStartTime,trackEndTime,costMatParam,...
        gapCloseParam,kalmanInfoLink,(1:numTracksLink)',...
        useLocalDensity.cg,nnDistLinkedFeat,useLocalDensity.nnWindowCG);

    %link tracks based on this cost matrix, allowing for birth and death
    [link12,link21] = lap(costMat,nonlinkMarker);
    link12 = double(link12);
    link21 = double(link21);

    %put the indices of all tracks from linking in one vector
    tracks2Link = (1:numTracksLink)';
    tracksRemaining = tracks2Link;
    
    %reserve memory space for matrix showing track connectivity
    compoundTrack = zeros(numTracksLink);

    %initialize compTrackIndx
    compTrackIndx = 0;

    while ~isempty(tracksRemaining)

        %update compound track index by 1
        compTrackIndx = compTrackIndx + 1;

        %take first track as a seed to build a compound track with
        %closed gaps and merges/splits
        trackSeed = tracksRemaining(1);
        seedLength = 1;
        seedLengthOld = 0; %dummy just to get into the while loop

        %while current seed contains more tracks than previous seed, i.e.
        %whie new track segments are still being added to the compound
        %track
        while seedLength > seedLengthOld

            %store current seed for later comparison
            seedLengthOld = seedLength;

            %find tracks connected to ends of seed tracks
            tmpTracks = link12(trackSeed);
            trackLink2End = tmpTracks(tmpTracks <= numTracksLink); %starts linked to ends
            trackMerge = [];
            if mergeSplit
                trackMerge = indxMerge(tmpTracks(tmpTracks > numTracksLink & ...
                    tmpTracks <= numTracksLink+numMerge) - numTracksLink); %tracks that ends merge with
            end

            %find tracks connected to starts of seed tracks
            tmpTracks = link21(trackSeed);
            trackLink2Start = tmpTracks(tmpTracks <= numTracksLink); %ends linked to starts
            trackSplit = [];
            if mergeSplit
                trackSplit = indxSplit(tmpTracks(tmpTracks > numTracksLink & ...
                    tmpTracks <= numTracksLink+numSplit) - numTracksLink); %tracks that starts split from
            end
            
            %put all tracks together as the new seed
            trackSeed = [trackSeed; trackLink2End; trackLink2Start; ...
                trackMerge; trackSplit];

            %remove repetitions and arrange tracks in ascending order
            trackSeed = unique(trackSeed);

            %get number of tracks in new seed
            seedLength = length(trackSeed);
            
            %expand new seed if merging/splitting are allowed
            if mergeSplit

                %variables storing merge/split seed tracks
                mergeSeed = [];
                splitSeed = [];
                
                %go over all seed tracks
                for iSeed = 1 : seedLength
                    
                    %get the location(s) of this track in indxMerge
                    mergeSeed = [mergeSeed; find(indxMerge == trackSeed(iSeed))];
                    
                    %get the location(s) of this track in indxSplit
                    splitSeed = [splitSeed; find(indxSplit == trackSeed(iSeed))];
                    
                end

                %add numTracksLink to mergeSeed and splitSeed to determine
                %their location in the cost matrix
                mergeSeed = mergeSeed + numTracksLink;
                splitSeed = splitSeed + numTracksLink;
                
                %find tracks merging with seed tracks
                trackMerge = [];
                for iSeed = 1 : length(mergeSeed)
                    trackMerge = [trackMerge; find(link12(1:numTracksLink)==mergeSeed(iSeed))];
                end

                %find tracks splitting from seed tracks
                trackSplit = [];
                for iSeed = 1 : length(splitSeed)
                    trackSplit = [trackSplit; find(link21(1:numTracksLink)==splitSeed(iSeed))];
                end

                %add these track to the seed
                trackSeed = [trackSeed; trackMerge; trackSplit];

                %remove repetitions and arrange tracks in ascending order
                trackSeed = unique(trackSeed);

                %get number of tracks in new seed
                seedLength = length(trackSeed);                
                
            end %(if mergeSplit)

        end %(while length(trackSeed) > length(trackSeedOld))

        %expand trackSeed to reserve memory for connetivity information
        trackSeedConnect = [trackSeed zeros(seedLength,2)];
        
        %store the tracks that the ends of the seed tracks are linked to,
        %and indicate whether it's an end-to-start link (+ve) or a merge (-ve)
        tmpTracks = link12(trackSeed);
        if mergeSplit
            tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                numTracksLink+numMerge) = -indxMerge(tmpTracks(tmpTracks > ...
                numTracksLink & tmpTracks <= numTracksLink+numMerge) - numTracksLink);
        end
        tmpTracks(tmpTracks > numTracksLink) = NaN;
        trackSeedConnect(:,2) = tmpTracks;
        
        %store the tracks that the starts of the seed tracks are linked to,
        %and indicate whether it's a start-to-end link (+ve) or a split (-ve)
        tmpTracks = link21(trackSeed);
        if mergeSplit
            tmpTracks(tmpTracks > numTracksLink & tmpTracks <= ...
                numTracksLink+numSplit) = -indxSplit(tmpTracks(tmpTracks > ...
                numTracksLink & tmpTracks <= numTracksLink+numSplit) - numTracksLink);
        end
        tmpTracks(tmpTracks > numTracksLink) = NaN;
        trackSeedConnect(:,3) = tmpTracks;

        %store tracks making up this compound track and their connectivity
        compoundTrack(compTrackIndx,1:3*seedLength) = reshape(...
            trackSeedConnect,3*seedLength,1)';

        %in the list of all tracks, indicate that these tracks have
        %been taken care of by placing NaN instead of their number
        tracks2Link(trackSeed) = NaN;

        %retain only tracks that have not been linked to anything yet
        tracksRemaining = tracks2Link(~isnan(tracks2Link));

    end %(while ~isempty(tracksRemaining))

    %remove empty rows
    maxValue = max(compoundTrack,[],2);
    compoundTrack = compoundTrack(maxValue > 0,:);

    %determine number of tracks after gap closing (including merge/split if
    %specified)
    numTracksCG = size(compoundTrack,1);

    %reserve memory for structure storing tracks after gap closing
    tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
        'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracksCG,1);

    %go over all compound tracks
    for iTrack = 1 : numTracksCG

        %get indices of tracks from linking making up current compound track
        %determine their number and connectivity
        trackSeedConnect = compoundTrack(iTrack,:)';
        trackSeedConnect = trackSeedConnect(trackSeedConnect ~= 0);
        seedLength = length(trackSeedConnect)/3; %number of segments making current track
        trackSeedConnect = reshape(trackSeedConnect,seedLength,3);
        
        %get their start times
        segmentStartTime = trackStartTime(trackSeedConnect(:,1));
        
        %arrange segments in ascending order of their start times
        [segmentStartTime,indxOrder] = sort(segmentStartTime);
        trackSeedConnect = trackSeedConnect(indxOrder,:);
        
        %get the segments' end times
        segmentEndTime = trackEndTime(trackSeedConnect(:,1));
        
        %calculate the segments' positions in the matrix of coordinates and
        %amplitudes
        segmentStartTime8 = 8 * (segmentStartTime - 1) + 1;
        segmentEndTime8   = 8 * segmentEndTime;

        %instead of having the connectivity in terms of the original track
        %indices, have it in terms of the indices of this subset of tracks
        %(which are arranged in ascending order of their start times)
        for iSeed = 1 : seedLength
            value = trackSeedConnect(iSeed,2);
            if value > 0
                trackSeedConnect(iSeed,2) = find(trackSeedConnect(:,1) == ...
                    value);
            elseif value < 0
                trackSeedConnect(iSeed,2) = -find(trackSeedConnect(:,1) == ...
                    -value);
            end
            value = trackSeedConnect(iSeed,3);
            if value > 0
                trackSeedConnect(iSeed,3) = find(trackSeedConnect(:,1) == ...
                    value);
            elseif value < 0
                trackSeedConnect(iSeed,3) = -find(trackSeedConnect(:,1) == ...
                    -value);
            end
        end

        %get track information from the matrices storing linking information
        tracksFeatIndxCG = tracksFeatIndxLink(trackSeedConnect(:,1),:);
        tracksCoordAmpCG = tracksCoordAmpLink(trackSeedConnect(:,1),:);

        %perform all gap closing links and modify connectivity accordingly
        %go over all starts in reverse order
        for iSeed = seedLength : -1 : 2
            
            %find the track this track might be connected to
            track2Append = trackSeedConnect(iSeed,3);

            %if there is a track (which is not a split)
            if track2Append > 0
                
                %put track information in the relevant row
                tracksFeatIndxCG(track2Append,segmentStartTime(iSeed):...
                    segmentEndTime(iSeed)) = tracksFeatIndxCG(iSeed,...
                    segmentStartTime(iSeed):segmentEndTime(iSeed));
                tracksFeatIndxCG(iSeed,:) = 0;
                tracksCoordAmpCG(track2Append,segmentStartTime8(iSeed):...
                    segmentEndTime8(iSeed)) = tracksCoordAmpCG(iSeed,...
                    segmentStartTime8(iSeed):segmentEndTime8(iSeed));
                tracksCoordAmpCG(iSeed,:) = NaN;
                
                %update segment information
                segmentEndTime(track2Append) = segmentEndTime(iSeed);
                segmentEndTime8(track2Append) = segmentEndTime8(iSeed);
                segmentEndTime(iSeed) = NaN;
                segmentEndTime8(iSeed) = NaN;
                segmentStartTime(iSeed) = NaN;
                segmentStartTime8(iSeed) = NaN;
                
                %update connectivity
                trackSeedConnect(track2Append,2) = trackSeedConnect(iSeed,2);
                trackSeedConnect(trackSeedConnect(:,2) == iSeed,2) = track2Append;
                trackSeedConnect(trackSeedConnect(:,3) == iSeed,3) = track2Append;
                trackSeedConnect(trackSeedConnect(:,2) == -iSeed,2) = -track2Append;
                trackSeedConnect(trackSeedConnect(:,3) == -iSeed,3) = -track2Append;

            end %(if track2Append > 0)
                
        end %(for iSeed = seedLength : -1 : 2)
        
        %find rows that are not empty
        maxValue = max(tracksFeatIndxCG,[],2);
        rowsNotEmpty = find(maxValue > 0);

        %remove empty rows
        tracksFeatIndxCG = tracksFeatIndxCG(rowsNotEmpty,:);
        tracksCoordAmpCG = tracksCoordAmpCG(rowsNotEmpty,:);
        segmentEndTime   = segmentEndTime(rowsNotEmpty);
        segmentStartTime = segmentStartTime(rowsNotEmpty);
        trackSeedConnect = trackSeedConnect(rowsNotEmpty,:);

        %update connectivity accordingly
        %by now, only merges and splits are left - thus no need for minus
        %sign to distinguish them from closed gaps
        for iSeed = 1 : length(rowsNotEmpty)
            trackSeedConnect(trackSeedConnect(:,2) == -rowsNotEmpty(...
                iSeed),2) = iSeed;
            trackSeedConnect(trackSeedConnect(:,3) == -rowsNotEmpty(...
                iSeed),3) = iSeed;
        end
        
        %determine new "seedLength"
        seedLength = length(rowsNotEmpty);

        %store the sequence of events of this track
        seqOfEvents = [segmentStartTime ones(seedLength,1) ...
            (1:seedLength)' trackSeedConnect(:,3); ...
            segmentEndTime 2*ones(seedLength,1) ...
            (1:seedLength)' trackSeedConnect(:,2)];
        
        %sort sequence of events in ascending order of time
        [tmp,indxOrder] = sort(seqOfEvents(:,1));
        seqOfEvents = seqOfEvents(indxOrder,:);
        
        %add 1 to the times of merges
        indx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2) == 2);
        seqOfEvents(indx,1) = seqOfEvents(indx,1) + 1;
        
        %store final tracks
        tracksFinal(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
        tracksFinal(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
        tracksFinal(iTrack).seqOfEvents = seqOfEvents;

    end %(for iTrack = 1 : numTracksCG)

end %(if any(trackStartTime > 1) && any(trackEndTime < numFrames))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([saveResDir filesep saveResFile],'costMatParam','gapCloseParam',...
    'kalmanInitParam','tracksFinal','tracksFeatIndxLink',...
    'tracksCoordAmpLink','kalmanInfoLink');


%%%%% ~~ the end ~~ %%%%%

