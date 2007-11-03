function [statsPerTrack] = scoreTrackLinksGapsLft(tracksFinal,tracksSim)
%SCORETRACKLINKSGAPSLFT evaluates links, gaps and lifetime per track
%
%SYNOPSIS [statsPerTrack] = scoreTrackLinksGapsLft(tracksFinal,tracksSim)
%
%INPUT  tracksFinal : Either output of trackCloseGapsKalman (structure) or
%                     output of trackWithGapClosing (matrix). Tracking must
%                     have been done on a simulated movieInfo (i.e. not
%                     obtained via detection).
%       tracksSim   : Simulated tracks as obtained from
%                     simulateMimickCD36_MS (structure format).
%
%OUTPUT statsPerTrack: (Number of tracks - 1) - by - 4 array. Row i
%                      corresponds to track i. The columns show:
%                      (1) Number of wrong links.
%                      (2) Number of wrongly closed gaps.
%                      (3) Lifetime of track.
%                      (4) Lifetime of corresponding ground truth track,
%                          when there is a ground truth track that starts
%                          at the same place and in the same frame. If
%                          there is no corresponding ground truth track,
%                          -1 will be shown instead.
%
%REMARKS I'm not sure how well the lifetime calculation works when there
%        are merges and splits.
%
%Khuloud Jaqaman, October 2007

%% process input variables

%extract track information out of tracksSim
xyCoordAll0 = trackInfoFromStruct(tracksSim);

if isstruct(tracksFinal)

    %extract track information out of tracksFinal
    xyCoordAll1 = trackInfoFromStruct(tracksFinal);

else

    %directly extract track information out of tracksFinal
    xyCoordAll1 = zeros(size(tracksFinal,1),size(tracksFinal,2)/4);
    xyCoordAll1(:,1:2:end) = tracksFinal(:,1:8:end);
    xyCoordAll1(:,2:2:end) = tracksFinal(:,2:8:end);
    xyCoordAll1(isnan(xyCoordAll1)) = 0;

end

%get "number of tracks" in tracking results
numTracks = size(xyCoordAll1,1);

%% compare links

%initialize output variable
statsPerTrack = NaN(numTracks,4);

%go over all tracks ...
for iTrack = 1 : numTracks

    %first, evaluate this track's links and closed gaps
    
    %initialize to zero number of wrong links and number of wrong gaps
    numLinksWrong = 0;
    numGapsWrong = 0;

    %find frames where this track exists
    framesExist = find(xyCoordAll1(iTrack,1:2:end));
    
    %go over these frames ...
    for iIndx = 1 : length(framesExist) - 1

        %get current frame and next frame
        iFrame1 = framesExist(iIndx);
        iFrame2 = framesExist(iIndx+1);

        %get coordinates in these two frames
        xyCoord12 = [xyCoordAll1(iTrack,2*iFrame1-1:2*iFrame1); ...
            xyCoordAll1(iTrack,2*iFrame2-1:2*iFrame2)];

        %find location of coordinates in first frame in ground truth
        [dummy,rowsGT] = intersect(xyCoordAll0(:,iFrame1*2-1:iFrame1*2),...
            xyCoord12(1,:),'rows');

        %get corresponding coordinates in second frame in ground truth
        xyCoord02 = [xyCoordAll0(rowsGT,iFrame2*2-1) xyCoordAll0(rowsGT,iFrame2*2)];

        %check whether they are the same as the coordinates in the second
        %frame in the tracking results
        commonCoord = intersect(xyCoord02,xyCoord12(2,:),'rows');

        %if the coordinates in the second frame are different between the
        %tracking results and the ground truth, then link/closed gap is wrong
        if isempty(commonCoord)
            if iFrame2 - iFrame1 == 1
                numLinksWrong = numLinksWrong + 1;
            else
                numGapsWrong = numGapsWrong + 1;
            end
        end

    end %(for iIndx = 1 : length(framesExist) - 1)
    
    %second, calculate the lifetime of this track
    trackLft = framesExist(end) - framesExist(1) + 1;
    
    %third, if the start of this track corresponds to the start of a track
    %in the ground truth, get the lifetime of the ground truth track
    
    %get coordinates in first frame of track
    iFrame1 = framesExist(1);
    xyCoord12 = xyCoordAll1(iTrack,2*iFrame1-1:2*iFrame1);

    %find location of coordinates in ground truth
    [dummy,rowsGT] = intersect(xyCoordAll0(:,iFrame1*2-1:iFrame1*2),...
        xyCoord12,'rows');
    
    %find frames where corresponding ground truth track(s) start(s)
    %and get their lifetimes if relevant
    groundTruthLft = -15 * ones(size(rowsGT));
    for iRow = 1 : length(rowsGT)
        framesExistGT = find(xyCoordAll0(rowsGT(iRow),1:2:end));
        if framesExistGT(1)==iFrame1
            groundTruthLft(iRow) = framesExistGT(end) - framesExistGT(1) + 1;
        end
    end
    
    %remove -1s which correspond to track(s) starting at a different time
    %from the track of interest
    groundTruthLft = groundTruthLft(groundTruthLft>0);
    
    if isempty(groundTruthLft) %if there is no corresponding track in GT, enter NaN
        groundTruthLft = -15;
    else %find GT track that has closest lifetime to track of interest
        tmp = abs(groundTruthLft - trackLft);
        groundTruthLft = groundTruthLft(tmp==min(tmp));
    end

    statsPerTrack(iTrack,:) = [numLinksWrong numGapsWrong trackLft ...
        groundTruthLft];

end %(for iTrack = 1 : numTracks)


%% Subfunction 1

function [xyCoordAll,xyCoordM,xyCoordS,trackStartRow] = trackInfoFromStruct(tracks)

%get number of tracks
numTracks = length(tracks);

%get number of frames
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%initialize output variables
xyCoordAll = [];
xyCoordM = [];
xyCoordS = [];
trackStartRow = zeros(numTracks,1);

%go over all tracks ...
for i = 1 : numTracks

    %get sequence of events of track
    seqOfEvents = tracks(i).seqOfEvents;

    %get start and end times of track
    startTime = seqOfEvents(1,1);
    endTime   = seqOfEvents(end,1);

    %extract track coordinates from structure
    tracksCoordAmpCG = tracks(i).tracksCoordAmpCG;
    xyCoordTmp = zeros(size(tracksCoordAmpCG,1),numFrames*2);
    xyCoordTmp(:,2*(startTime-1)+1:2:2*(endTime-1)+1) = tracksCoordAmpCG(:,1:8:end);
    xyCoordTmp(:,2*(startTime-1)+2:2:2*(endTime-1)+2) = tracksCoordAmpCG(:,2:8:end);

    %determine row where this compound track starts in xyCoordAll
    trackStartRow(i) = size(xyCoordAll,1) + 1;

    %get merge events
    mergeIndx = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)));

    %go over all merges
    for iMerge = mergeIndx'

        %get merge time
        timeMerge = seqOfEvents(iMerge,1);

        %store coordinates belonging to this merge in xyCoordM
        xyCoordM = [xyCoordM; [timeMerge ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-1)+1:2*(timeMerge-1)+2) ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-2)+1:2*(timeMerge-2)+2) ...
            xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-2)+1:2*(timeMerge-2)+2)]];

        %add coordinates after merge to merging track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iMerge,3),2*(timeMerge-1)+1:2*(timeMerge-1)+2) = ...
            xyCoordTmp(seqOfEvents(iMerge,4),2*(timeMerge-1)+1:2*(timeMerge-1)+2);

    end

    %get split events
    splitIndx = find(seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)));

    %go over all splits
    for iSplit = splitIndx'

        %get split time
        timeSplit = seqOfEvents(iSplit,1);

        %store coordinates belonging to this split in xyCoordS
        xyCoordS = [xyCoordS; [timeSplit ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2) ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-1)+1:2*(timeSplit-1)+2) ...
            xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-1)+1:2*(timeSplit-1)+2)]];

        %add coordinates before split to splitting track in xyCoordTmp
        xyCoordTmp(seqOfEvents(iSplit,3),2*(timeSplit-2)+1:2*(timeSplit-2)+2) = ...
            xyCoordTmp(seqOfEvents(iSplit,4),2*(timeSplit-2)+1:2*(timeSplit-2)+2);

    end

    %store coordinates in xyCoordAll
    xyCoordAll = [xyCoordAll; xyCoordTmp];

end

%replace NaNs by zero
xyCoordAll(isnan(xyCoordAll)) = 0;
