function [mergesInfo,splitsInfo,mergesInfoSpace,splitsInfoSpace] = ...
    findMergesSplits(tracks,probDim,removePotArtifacts,plotRes,calcTrackType)
%FINDMERGESSPLITS finds the merges and splits in each track and gives back their time and location
%
%SYNOPSIS [mergesInfo,splitsInfo,mergesInfoSpace,splitsInfoSpace] = ...
%    findMergesSplits(tracks,probDim,removePotArtifacts,plotRes)
%
%INPUT  tracks    : Output of trackCloseGapsKalman:
%                   Structure array with number of entries equal to
%                   the number of tracks (or compound tracks when
%                   merging/splitting are considered). Contains the
%                   fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       probDim   : 2 for 2D, 3 for 3D. Optional. Default: 2.
%       removePotArtifacts: 1 to remove potentially artifactual merges and
%                   splits, resulting for instance from detection
%                   artifacts, 0 otherwise. 
%                   Optional. Default: 1.
%       plotRes   : 0 to not plot anything, 1 to make a spatial map of
%                   merges and splits.
%                   Optional. Default: 0.
%       calcTrackType: Estimate track type based on asymmetry. 1 to
%                   estimate, 0 to not estimate. 
%                   Optional. Default: 1.
%
%OUTPUT mergesInfo     : 2D array where first column indicates track number,
%                        second column indicates track type (1 linear, 0 o.w.),
%                        third column indicates number of merges, and 
%                        subsequent columns indicate merge times.
%                        Track type will be 0 if calcTrackType = 0.
%       splitsInfo     : 2D array where first column indicates track number,
%                        second column indicates track type (1 linear, 0 o.w.),
%                        third column indicates number of splits, and 
%                        subsequence columns indicate split times.
%                        Track type will be 0 if calcTrackType = 0.
%       mergesInfoSpace: 2D array that is a continuation of mergesInfoTime,
%                        storing the (x,y,[z])-coordinates of each merge.
%                        Every row corresponds to the same row in
%                        mergesInfo. Every merge gets 2 (in 2D) or 3 (in
%                        3D) columns for x, y and z (if 3D).
%       splitsInfoSpace: 2D array that is a continuation of splitsInfoTime,
%                        storing the (x,y,[z])-coordinates of each split.
%                        Every row corresponds to the same row in
%                        splitsInfo. Every split gets 2 (in 2D) or 3 (in
%                        3D) columns for x, y and z (if 3D).
%
%
%REMARKS Plotting implemented for 2D only.
%
%Khuloud Jaqaman, October 2007

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--findMergesSplits: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(probDim)
    probDim = 2;
end

if nargin < 3 || isempty(removePotArtifacts)
    removePotArtifacts = 1;
end

if nargin < 4 || isempty(plotRes)
    plotRes = 0;
end
if probDim ~= 2
    plotRes = 0;
end

if nargin < 5 || isempty(calcTrackType)
    calcTrackType = 1;
end

%get number of tracks
numTracks = length(tracks);

%get number of segments per track
numSegments = getNumSegmentsPerTrack(tracks);

%estimate track types
if calcTrackType
    trackType = getTrackType(tracks,probDim);
else
    trackType = zeros(numTracks,1);
end

% %convert tracks from structure to matrix format
% tracksMat = convStruct2MatIgnoreMS(tracks);
% 
% %find minimum and maximum coordinates for plotting
% xCoord = tracksMat(:,1:8:end);
% minXCoord = min(floor(min(xCoord(:))),0);
% maxXCoord =  ceil(max(xCoord(:)));
% yCoord = tracksMat(:,2:8:end);
% minYCoord = min(floor(min(yCoord(:))),0);
% maxYCoord =  ceil(max(yCoord(:)));

%% Merge/split statistics

[mergesInfo,splitsInfo] = deal(zeros(numTracks,max(numSegments)));
[mergesInfoSpace,splitsInfoSpace] = deal(repmat(mergesInfo,1,3));

%go over all tracks ...
for iTrack = 1 : numTracks
    
    %get track's sequence of events
    seqOfEvents = tracks(iTrack).seqOfEvents;
    
    %remove splits and merges that are most likely artifacts
    if removePotArtifacts
        seqOfEvents = removeSplitMergeArtifacts(seqOfEvents,0);
    end

    %get track's coordinates
    trackCoordX = tracks(iTrack).tracksCoordAmpCG(:,1:8:end);
    trackCoordY = tracks(iTrack).tracksCoordAmpCG(:,2:8:end);
    trackCoordZ = tracks(iTrack).tracksCoordAmpCG(:,3:8:end);
    
    %find rows with merging information
    indxMerge = find( seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)) );
    
    %get the merge times (absolute and relative to track start)
    mergeTimes = seqOfEvents(indxMerge,1);
    relMergeTimes = mergeTimes - seqOfEvents(1,1) + 1;
    
    %get the track segments that are merged with
    mergeSegment = seqOfEvents(indxMerge,4);
    
    %get the coordinates of each merge
    mergeCoords = [];
    for iMerge = 1 : length(indxMerge)
        mergeCoords = [mergeCoords ...
            trackCoordX(mergeSegment(iMerge),relMergeTimes(iMerge)) ...
            trackCoordY(mergeSegment(iMerge),relMergeTimes(iMerge)) ...
            trackCoordZ(mergeSegment(iMerge),relMergeTimes(iMerge))]; %#ok<AGROW>
    end
    
    %store the merge information for this track
    mergesInfo(iTrack,1:length(mergeTimes)+2) = [trackType(iTrack) ...
        length(mergeTimes) mergeTimes'];
    mergesInfoSpace(iTrack,1:3*length(mergeTimes)) = mergeCoords;
        
    %find rows with splitting information
    indxSplit = find( seqOfEvents(:,2)==1 & ~isnan(seqOfEvents(:,4)) );
    
    %get split times (absolute and relative to track start)
    splitTimes = seqOfEvents(indxSplit,1);
    relSplitTimes = splitTimes - seqOfEvents(1,1) + 1;

    %get the track segments that are split from
    splitSegment = seqOfEvents(indxSplit,4);
    
    %get the coordinates of each split
    splitCoords = [];
    for iSplit = 1 : length(indxSplit)
        splitCoords = [splitCoords ...
            trackCoordX(splitSegment(iSplit),relSplitTimes(iSplit)-1) ...
            trackCoordY(splitSegment(iSplit),relSplitTimes(iSplit)-1) ...
            trackCoordZ(splitSegment(iSplit),relSplitTimes(iSplit)-1)]; %#ok<AGROW>
    end
    
    %store the split information for this track
    splitsInfo(iTrack,1:length(splitTimes)+2) = [trackType(iTrack) ...
        length(splitTimes) splitTimes'];
    splitsInfoSpace(iTrack,1:3*length(splitTimes)) = splitCoords;

end

%remove empty columns
fullIndx = find(sum(mergesInfo(:,2:end))~=0);
mergesInfo = [mergesInfo(:,1) mergesInfo(:,1+fullIndx)];
fullIndx = sum(mergesInfoSpace)~=0;
mergesInfoSpace = mergesInfoSpace(:,fullIndx);
fullIndx = find(sum(splitsInfo(:,2:end))~=0);
splitsInfo = [splitsInfo(:,1) splitsInfo(:,1+fullIndx)];
fullIndx = sum(splitsInfoSpace)~=0;
splitsInfoSpace = splitsInfoSpace(:,fullIndx);

%remove rows without merges or splits
filledRows = find(any(mergesInfo(:,2)~=0,2));
mergesInfo = [filledRows mergesInfo(filledRows,:)];
mergesInfoSpace = mergesInfoSpace(filledRows,:);
filledRows = find(any(splitsInfo(:,2)~=0,2));
splitsInfo = [filledRows splitsInfo(filledRows,:)];
splitsInfoSpace = splitsInfoSpace(filledRows,:);

%% Plotting

if plotRes
    switch probDim
        case 2
            plotMergeSplitPositions2D(tracks,mergesInfo,splitsInfo,...
                mergesInfoSpace,splitsInfoSpace)
    end
end

%% %%%%% ~~ the end ~~ %%%%%

