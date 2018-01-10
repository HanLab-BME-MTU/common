function [probMotionMode,modeMotionChar,errFlag] = summarizeDiffModeRes(tracks,...
    diffAnalysisRes,numMode,minTrackLen,probDim,extractType)
%SUMMARIZEDIFFMODERES calculates diffusion mode probabilities and motion characteristics from diffusion mode analysis
%
%SYNOPSIS [probMotionMode,modeMotionChar,errFlag] = summarizeDiffModeRes(tracks,...
%    minTrackLen,probDim,diffAnalysisRes,extractType,numMode)    
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       diffAnalysisRes: Diffusion mode analysis results (output of
%                    trackDiffModeAnalysis).
%       numMode    : Number of modes.
%       minTrackLen: Minimum length of a track to be used in getting
%                    motion type statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%       extractType: 1 - Analyze every track segment separately.
%                    2 - Extract from each compound track the longest
%                        trajectory to use in analysis - NOT IMPLEMENTED
%                        YET.
%                    Must use same extractType as in trackDiffModeAnalysis.
%                    Variable irrelevant if tracks are input as a matrix.
%                    Optional. Default: 1.
%
%OUTPUT probMotionMode: (Number of modes + 1) - by - 2 array. 
%                    Rows refer to mode 1, mode 2, mode 3, etc.
%                    Last row is for tracks with undetermined mode (too
%                    short for analysis).
%                    1st column: Each mode's probability.
%                    2rd column: Number of features per frame in each mode.
%       modeMotionChar:Structure array summarizing motion characteristics
%                    for each diffusion mode. Last entry is for
%                    unclassified tracks.
%                    It contains three fields:
%           .distribution: Nx3 array. Each row belongs to a separate
%                          trajectory in the category. The columns
%                          store the diffusion coefficient,
%                          localization precision, and trajectory
%                          lifetime.
%           .meanStd     : 2x3 array. Columns are same as "distribution".
%                          1st row is the mean, 2nd row is the standard
%                          deviation.
%           .median      : 1x3 array showing the median of the columns
%                          in "distribution".
%
%Khuloud Jaqaman, November 2017

%% output

probMotionMode = [];
modeMotionChar = [];
errFlag = 0;

%% input

if nargin < 3 || isempty(tracks) || isempty(diffAnalysisRes) || isempty(numMode)
    disp('summarizeDiffModeRes: Missing input arguments!');
    return
end

if nargin < 4 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 5 || isempty(probDim)
    probDim = 2;
end

if nargin < 6 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--summarizeDiffModeRes: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
if ~isempty(indx)
    tracks = tracks(indx);
    diffAnalysisRes = diffAnalysisRes(indx);
else
    disp('Ignoring minTrackLen because imposing it leaves no tracks to analyze.')
end

%get number of tracks and number of frames
numTracks = length(tracks);
seqOfEvents = vertcat(tracks.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%put tracks in matrix format
if extractType == 1
    [tracksMat,tracksIndxMat,trackStartRow] = convStruct2MatIgnoreMS(tracks);
end

%get number of track segments
numTrackSegments = size(tracksMat,1);

%get track lengths
trackSegmentLft = getTrackSEL(tracksMat);
trackSegmentLft = trackSegmentLft(:,3);

%% features

%get average number of features per frame
numFeatTot = length(find(tracksIndxMat(:)));
aveFeatPerFrame = numFeatTot / numFrames;

%% motion mode probabilities

%get track segment modes from diffusion mode analysis
trackSegmentMode = vertcat(diffAnalysisRes.diffMode);
numModes = numMode + 1; %add 1 to account for unclassified tracks

%get indices of the different modes
indxMode = cell(numModes,1);
for i = 1 : numModes - 1
    indxMode{i} = find(trackSegmentMode == i);
end
indxMode{end} = find(isnan(trackSegmentMode)); %unclassified

%calculate number of track segments per mode
%calculate number of features in each motion mode
[numSegmentsMode,numFeatMode] = deal(zeros(numModes,1));
for i = 1 : numModes
    numSegmentsMode(i) = length(indxMode{i});
    numFeatMode(i) = length(find(tracksIndxMat(indxMode{i},:)));
end

%calculate fraction of track segments falling in each mode
fracSegmentsMode = numSegmentsMode / numTrackSegments;

%get fraction of features undergoing each motion mode - this is the
%probability of a feature undergoing a certain motion mode
probFeatMode = numFeatMode / numFeatTot;

%put output together
probMotionMode = [probFeatMode numFeatMode/numFrames];

%% motion characteristics

%initialize structure
modeMotionChar = repmat(struct('distribution',[],'meanStd',[],'median',[]),numModes,1);

%extract diffusion coefficients and positional standard deviation
diffCoefAll = vertcat(diffAnalysisRes.diffCoef);
meanPosStdAll = vertcat(diffAnalysisRes.meanPosStd);

%distribute motion characteristics based on motion mode
for iMode = 1 : numModes
    motionCharTmp = [diffCoefAll(indxMode{iMode}) meanPosStdAll(indxMode{iMode}) trackSegmentLft(indxMode{iMode})];
    modeMotionChar(iMode).distribution = motionCharTmp;
    modeMotionChar(iMode).meanStd = [nanmean(motionCharTmp,1); nanstd(motionCharTmp,[],1)];
    modeMotionChar(iMode).median = nanmedian(motionCharTmp,1);
end

%% ~~~ the end ~~~

