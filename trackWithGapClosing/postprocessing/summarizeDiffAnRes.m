function [probMotionType,motionChar,errFlag] = summarizeDiffAnRes(tracks,...
    minTrackLen,probDim,diffAnalysisRes,extractType)
%SUMMARIZEDIFFANRES calculates motion type probabilities and motion characteristics from diffusion analysis
%
%SYNOPSIS [probMotionType,motionChar,errFlag] = summarizeDiffAnRes(tracks,...
%    minTrackLen,probDim,diffAnalysisRes,extractType)
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       minTrackLen: Minimum length of a track to be used in getting
%                    motion type statistics.
%                    Optional. Default: 5.
%       probDim    : Dimensionality - 2 for 2D, 3 for 3D.
%                    Optional. Default: 2.
%       diffAnalysisRes: Diffusion analysis results (output of
%                    trackDiffAnalysis1). Optional. If not input, it will
%                    be calculated.
%       extractType: 1 - Analyze every track segment separately.
%                    2 - Extract from each compound track the longest
%                        trajectory to use in analysis - NOT IMPLEMENTED
%                        YET.
%                    Must use same extractType as in trackDiffusionAnalysis1.
%                    Variable irrelevant if tracks are input as a matrix.
%                    Optional. Default: 1.
%
%OUTPUT probMotionType: 9-by-2 array. Rows refer to:
%                    (1) Linear & 1D confined,
%                    (2) Linear & 1D Brownian,
%                    (3) Linear & 1D directed,
%                    (4) Linear & diffusion undetermined,
%                    (5) Not linear/unclassied & 2D confined,
%                    (6) Not linear/unclassified & 2D Brownian,
%                    (7) Not linear/unclassified & 2D directed,
%                    (8) Not linear & diffusion undetermined, and
%                    (9) Completely undetermined.
%                    1st column: Each motion type's probability.
%                    2nd column: Probability of a motion type within its
%                    category, i.e. 1-4 are in the linear category, while
%                    1-8 or 1-7&9 are in the not-linear category (1-8 are
%                    used when asymmetry is checked for, 1-7&9 are used
%                    when asymmetry is not checked for).
%
%Khuloud Jaqaman, March 2010

%% output

probMotionType = [];
motionChar = [];
errFlag = 0;

%% input

if nargin < 1 || isempty(tracks)
    disp('summarizeDiffAnRes: Missing input argument!');
    return
end

if nargin < 2 || isempty(minTrackLen)
    minTrackLen = 5;
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

if nargin < 4 || isempty(diffAnalysisRes)
    [diffAnalysisRes,errFlag] = trackDiffusionAnalysis1(tracks,1,probDim,...
        1,[0.05 0.05],0);
    if errFlag
        return
    end
end

if nargin < 5 || isempty(extractType)
    extractType = 1;
else
    if ~any(extractType == [1 2])
        disp('--trackDiffusionAnalysis1: Variable extractType should be 1 or 2.');
        errFlag = 1;
    end
end

%% preamble

%keep only tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
indx = chooseTracks(tracks,criteria);
clear criteria
tracks = tracks(indx);
diffAnalysisRes = diffAnalysisRes(indx);

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

%% motion type probabilities

%get track segment classification from diffusion analysis
trackSegmentClass = vertcat(diffAnalysisRes.classification);

%get indices of the different motion types
indxLinConf     = find( trackSegmentClass(:,1) == 1   & trackSegmentClass(:,3) == 1   );
indxLinBrown    = find( trackSegmentClass(:,1) == 1   & trackSegmentClass(:,3) == 2   );
indxLinDir      = find( trackSegmentClass(:,1) == 1   & trackSegmentClass(:,3) == 3   );
indxLinUndet    = find( trackSegmentClass(:,1) == 1   & isnan(trackSegmentClass(:,3)) );
indxNonlinConf  = find( trackSegmentClass(:,1) ~= 1   & trackSegmentClass(:,2) == 1   );
indxNonlinBrown = find( trackSegmentClass(:,1) ~= 1   & trackSegmentClass(:,2) == 2   );
indxNonlinDir   = find( trackSegmentClass(:,1) ~= 1   & trackSegmentClass(:,2) == 3   );
indxNonlinUndet = find( trackSegmentClass(:,1) == 0   & isnan(trackSegmentClass(:,2)) );
indxUndetUndet  = find( isnan(trackSegmentClass(:,1)) & isnan(trackSegmentClass(:,2)) );

%calculate number of track segments per motion type
numSegmentsType = [length(indxLinConf) length(indxLinBrown) ...
    length(indxLinDir) length(indxLinUndet) ...
    length(indxNonlinConf) length(indxNonlinBrown) ...
    length(indxNonlinDir) length(indxNonlinUndet) length(indxUndetUndet)]';

%calculate fraction of track segments falling in each motion type
fracSegmentsType = numSegmentsType / numTrackSegments;

%calculate number of features in each motion type
numFeatType = [length(find(tracksIndxMat(indxLinConf,:))) ...
    length(find(tracksIndxMat(indxLinBrown,:))) ...
    length(find(tracksIndxMat(indxLinDir,:))) ...
    length(find(tracksIndxMat(indxLinUndet,:))) ...
    length(find(tracksIndxMat(indxNonlinConf,:))) ...
    length(find(tracksIndxMat(indxNonlinBrown,:))) ...
    length(find(tracksIndxMat(indxNonlinDir,:))) ...
    length(find(tracksIndxMat(indxNonlinUndet,:))) ...
    length(find(tracksIndxMat(indxUndetUndet,:)))]';

%get fraction of features undergoing each motion type - this is the
%probability of a feature undergoing a certain motion type
probFeatType = numFeatType / numFeatTot;

%within the overall linear and non-linear motion types, calculate the
%probabilities of the different sub-types
probFeatLinSubType = probFeatType(1:4) / sum(probFeatType(1:4));
if all(isnan(trackSegmentClass(:,1)))
    probFeatNonlinSubType = probFeatType([5:7 9]) / sum(probFeatType([5:7 9]));
    probFeatSubType = [probFeatLinSubType; probFeatNonlinSubType(1:3); NaN; ...
        probFeatNonlinSubType(4)];
else
    probFeatNonlinSubType = probFeatType(5:8) / sum(probFeatType(5:8));
    probFeatSubType = [probFeatLinSubType; probFeatNonlinSubType; NaN];
end

%combine the motion type and sub-type probabilities
probMotionType = [probFeatType probFeatSubType];

%% motion characteristics

%extract diffusion coefficients and confinement radii
diffCoefAll = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');
confRadAll = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius');

%distribute motion characteristics based on motion types

%linear, all together
motionCharTmp = [diffCoefAll([indxLinConf;indxLinBrown;indxLinDir]) ...
    confRadAll([indxLinConf;indxLinBrown;indxLinDir],:) ...
    trackSegmentLft([indxLinConf;indxLinBrown;indxLinDir])];
motionChar.linear.all.distribution = motionCharTmp;
motionChar.linear.all.meanStd = [nanmean(motionCharTmp); nanstd(motionCharTmp)];

%non-linear & confined
motionCharTmp = [diffCoefAll(indxNonlinConf) confRadAll(indxNonlinConf,:) ...
    trackSegmentLft(indxNonlinConf)];
motionChar.notLinear.confined.distribution = motionCharTmp;
motionChar.notLinear.confined.meanStd = [nanmean(motionCharTmp); nanstd(motionCharTmp)];

%non-linear & Brownian
motionCharTmp = [diffCoefAll(indxNonlinBrown) confRadAll(indxNonlinBrown,:) ...
    trackSegmentLft(indxNonlinBrown)];
motionChar.notLinear.brownian.distribution = motionCharTmp;
motionChar.notLinear.brownian.meanStd = [nanmean(motionCharTmp); nanstd(motionCharTmp)];

%non-linear & directed
motionCharTmp = [diffCoefAll(indxNonlinDir) confRadAll(indxNonlinDir,:) ...
    trackSegmentLft(indxNonlinDir)];
motionChar.notLinear.directed.distribution = motionCharTmp;
motionChar.notLinear.directed.meanStd = [nanmean(motionCharTmp); nanstd(motionCharTmp)];

%% ~~~ the end ~~~

