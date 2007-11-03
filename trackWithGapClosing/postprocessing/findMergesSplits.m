function [mergesInfo,splitsInfo] = findMergesSplits(tracks)
%FINDMERGESSPLITS finds the merges and splits in each track and gives back their location
%
%SYNOPSIS [mergesInfo,splitsInfo] = findMergesSplits(tracks)
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
%
%OUTPUT mergesInfo: 2D array where first column indicates number of merges
%                   per track and subsequent columns indicate merge times.
%       splitsInfo: 2D array where first column indicates number of splits
%                   per track and subsequence columns indicate split times.
%
%Khuloud Jaqaman, October 2007

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--tracks: Incorrect number of input arguments!');
    return
end

%get number of tracks
numTracks = length(tracks);

%get number of segments per track
numSegments = zeros(numTracks,1);
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracks(iTrack).tracksCoordAmpCG,1);
end

%% Merge/split stats

mergesInfo = zeros(numTracks,max(numSegments));
splitsInfo = mergesInfo;

%go over all tracks ...
for iTrack = 1 : numTracks
    
    %get track's sequence of events
    seqOfEvents = tracks(iTrack).seqOfEvents;
    
    %get number of merges for this track
    mergeTimes = seqOfEvents((seqOfEvents(:,2)==2&~isnan(seqOfEvents(:,4))),1);
    mergesInfo(iTrack,1:length(mergeTimes)+1) = [length(mergeTimes) mergeTimes'];
    
    %get number of splits for this track
    splitTimes = seqOfEvents((seqOfEvents(:,2)==1&~isnan(seqOfEvents(:,4))),1);
    splitsInfo(iTrack,1:length(splitTimes)+1) = [length(splitTimes) splitTimes'];

end

%remove empty columns
mergesInfo = mergesInfo(:,any(mergesInfo~=0,1));
splitsInfo = splitsInfo(:,any(splitsInfo~=0,1));

%remove empty rows
filledRows = find(any(mergesInfo~=0,2));
mergesInfo = [filledRows mergesInfo(filledRows,:)];
filledRows = find(any(splitsInfo~=0,2));
splitsInfo = [filledRows splitsInfo(filledRows,:)];


%% %%%%% ~~ the end ~~ %%%%%

