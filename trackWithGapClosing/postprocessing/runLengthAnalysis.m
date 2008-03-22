function [runLengthPos,runLengthNeg,dispPos,dispNeg,positionProj] = ...
    runLengthAnalysis(tracks,centerCoord,minTrackLen)
%MSSTATS calculate some merge/split statistics
%
%SYNOPSIS [runLengthPos,runLengthNeg,dispPos,dispNeg,positionProj] = ...
%    runLengthAnalysis(tracks,centerCoord,minTrackLen)
%
%INPUT  tracks     : Output of trackCloseGapsKalman.
%       centerCoord: Coordinates of convergence point of tracks ("center
%                    of nucleus").
%       minTrackLen: Minimum length of a track to be used in analysis.
%                    Optional. Default: 5.
%
%OUTPUT runLengthPos: Distribution of run lengths (i.e. number of
%                     consecutive steps in the same direction) when moving
%                     away from the nucleus.
%       runLengthNeg: Distribution of run lengths when moving toward the
%                     nucleus.
%       dispPos     : Distribution of displacement magntidues when moving
%                     away from the nucleus.
%       dispNeg     : Distribution of displacement magnitudes when moving
%                     toward the nucleus.
%       positionProj: Structure array of position projection onto 
%                     preferred direction of motion, for input into
%                     armaxFitKalman (can be converted into input for
%                     trajectoryAnalysis using convertTrajectoryData).
%
%Khuloud Jaqaman, February 2008

%% output
runLengthPos = [];
runLengthNeg = [];
dispPos = [];
dispNeg = [];

%% input

if nargin < 2 || isempty(tracks) || isempty (centerCoord)
    disp('calcStatsMS: Missing input arguments!');
    return
end

if nargin < 3 || isempty(minTrackLen)
    minTrackLen = 5;
end

%% preamble

%ignore merges and splits and divide compound tracks back into the
%individual tracks
inputStruct = tracks;
clear tracksFinal
tracks = convStruct2MatIgnoreMS(inputStruct);

%keep only linear tracks with length >= minTrackLen
criteria.lifeTime.min = minTrackLen;
criteria.trackType = 1;
indx = chooseTracks(tracks,criteria);
tracks = tracks(indx,:);
numTracks = size(tracks,1);

%reserve memory for displacement projection structure
positionProj = repmat(struct('observations',[],'time',[]),numTracks,1);

%% displacement statistics

%go over all tracks ...
for iTrack = 1 : numTracks

    %get the positions in this track and their standard deviations
    %keep NaNs to mark gaps
    trackCoordX = tracks(iTrack,1:8:end)';
    deltaCoordX = tracks(iTrack,5:8:end)';
    trackCoordY = tracks(iTrack,2:8:end)';
    deltaCoordY = tracks(iTrack,6:8:end)';
    trackCoord = [trackCoordX trackCoordY];
    deltaCoord = [deltaCoordX deltaCoordY];
    
    %project positions onto track's direction of motion
    [posAlongDir,deltaPosAlongDir,velDir] = projectCoordOntoDir(...
        trackCoord,deltaCoord,[],centerCoord);

    %store track information in dispProj
    positionProj(iTrack).observations = [posAlongDir deltaPosAlongDir];
    positionProj(iTrack).time = [(1:length(posAlongDir))' zeros(length(posAlongDir),1)];
    
    %calculate vector of displacements
    %keep NaNs to mark gaps
    dispVec = trackCoord(2:end,:) - trackCoord(1:end-1,:);
    
    %calculate the dot product of displacements with the direction vector
    dispAlongDir = dispVec * velDir;

    %separate positive displacements from negative displacements
    %NaNs (gaps) have no effect here
    dispPosT = dispAlongDir(dispAlongDir>0);
    dispNegT = -dispAlongDir(dispAlongDir<0);
    
    %get the sign of displacements (regardless of value)
    %gaps will constribute NaNs here
    dispSign = sign(dispAlongDir);
    
    %take the difference to find transition points
    dispSignDiff = diff(dispSign);
    transValue = [0; dispSignDiff(dispSignDiff~=0)];
    transPoints = [0; find(dispSignDiff)];
    
    %determine run length, i.e. consecutive steps in each direction before
    %switching
    runLength = diff(transPoints);
    
    %store run lengths away from nucleus (pos) separate from run lengths
    %toward nucleus (neg)
    %remove run lengths involving gaps (NaNs)
    runLengthPosT = [];
    runLengthNegT = [];
    for iTrans = 1 : length(runLength)
        if ~isnan(transValue(iTrans+1)) && ~isnan(transValue(iTrans))
            if transValue(iTrans+1) > 0
                runLengthNegT = [runLengthNegT; runLength(iTrans)];
            else
                runLengthPosT = [runLengthPosT; runLength(iTrans)];
            end
        end
    end
    
    %separate positive steps from negative steps
    %add this track's results to the rest
    runLengthPos = [runLengthPos; runLengthPosT];
    runLengthNeg = [runLengthNeg; runLengthNegT];
    dispPos = [dispPos; dispPosT];
    dispNeg = [dispNeg; dispNegT];
    
end

%% ~~~ the end ~~~



