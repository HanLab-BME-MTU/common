function [msStats,mergesInfo,splitsInfo] = calcStatsMS(tracks,minTrackLen,probDim)
%MSSTATS calculate some merge/split statistics
%
%SYNOPSIS [msStats,mergesInfo,splitsInfo] = calcStatsMS(tracks,minTrackLen,probDim)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    merge/split statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%
%OUTPUT msStats    : Row vector with entries: 
%                    1st: number of features per frame; 
%                    2nd/3rd: number of merges/splits per feature;
%                    4th: total number of tracks (length >= minTrackLen);
%                    5th: fraction of tracks that are linear;
%                    6th: fraction of features undergoing linear motion;
%                    7th/8th: fraction of merges/splits in linear tracks.
%       mergesInfo : Output of findMergesSplits;
%       splitsInfo : Output of findMergesSplits;
%
%Khuloud Jaqaman, December 2007

%% input

if nargin < 1 || isempty(tracks)
    disp('msStats: Missing input argument!');
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%% features

%get average number of features per frame
numFeatTot = 0; 
for iTrack = 1 : numTracks
    numFeatTot = numFeatTot + length(find(tracks(iTrack).tracksFeatIndxCG)); 
end
aveFeatPerFrame = numFeatTot / numFrames;

%% track types

%find number of tracks per type
criteria.trackType = 0;
indxBrown = chooseTracks(tracks,criteria);
numTracksBrown = length(indxBrown);
criteria.trackType = 1;
indxLin = chooseTracks(tracks,criteria);
numTracksLin = length(indxLin);

%calculate fraction of tracks that are linear
fracTrackLin = numTracksLin / (numTracksBrown + numTracksLin);

%calculate number of features undergoing linear motion
numFeatLin = 0;
for iTrack = indxLin'
    numFeatLin = numFeatLin + length(find(tracks(iTrack).tracksFeatIndxCG));
end

%get fraction of features undergoing linear motion
fracFeatLin = numFeatLin / numFeatTot;

%% merges/splits vs. features

%locate merges and splits in tracks
[mergesInfo,splitsInfo] = findMergesSplits(tracks,probDim);

%get total number of merges and splits
numMergesTot = sum(mergesInfo(:,3));
numSplitsTot = sum(splitsInfo(:,3));

%calculate average number of merges/splits per frame
aveMergePerFrame = numMergesTot / numFrames;
aveSplitPerFrame = numSplitsTot / numFrames;

%calculate average number of merges/splits per feature
aveMergePerFeat = aveMergePerFrame / aveFeatPerFrame;
aveSplitPerFeat = aveSplitPerFrame / aveFeatPerFrame;

%% merges/splits vs. track type

%get number of merges/splits happening in linear tracks
numMergesLin = sum(mergesInfo(mergesInfo(:,2)==1,3));
numSplitsLin = sum(splitsInfo(splitsInfo(:,2)==1,3));

%get their fraction of total number
fracMergesLin = numMergesLin / numMergesTot;
fracSplitsLin = numSplitsLin / numSplitsTot;

%% output

msStats = [aveFeatPerFrame aveMergePerFeat aveSplitPerFeat ...
    numTracks fracTrackLin fracFeatLin fracMergesLin fracSplitsLin];

%% ~~~ the end ~~~



