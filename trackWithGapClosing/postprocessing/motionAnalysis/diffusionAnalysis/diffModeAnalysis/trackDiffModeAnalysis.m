function diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStruct)
%TRACKDIFFMODEANALYSIS classifies tracks into diffusion modes
%
%SYNOPSIS diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStruct)
%
%INPUT  tracksFinal : Output of trackCloseGapsKalman.
%       diffModeDividerStruct: Structure with fields:
%           .trajLength: Trajectory lengths for which mode divider values
%                        are supplied.
%           .coordStd  : Coordinate standard deviations for which mode
%                        divider values are supplied.
%           .divider   : (number of trajectory lengths) x (number of
%                        coordinates stds) x (number of modes - 1) array
%                        of mode divider values.
%                     Trajectories shorter than the minimum length with a
%                     divider cannot be classified.
%                     Trajectories longer than the maximum length with a
%                     divider will be classified using the divider values
%                     for the maximum length.
%                     The coordStd used is the closest to each trajectory's
%                     average coordinate standard deviation.
%                     *** If no struct or [] is input, then trajectories are not
%                     classified into modes ***
%
%OUTPUT diffModeAnalysisRes : Structure array with the following fields per
%                         track:
%           .diffMode: Diffusion mode.
%           .diffCoef: Diffusion coefficient as calculated from the mean
%           square frame-to-frame displacement and mean positional standard
%           deviaton.
%           .msdF2F  : Mean square frame-to-frame displacement.
%           .meanPosStd: Mean positonal standard deviation.
%
%REMARKS Code is written for 2D case only, but can be generalized to 3D.
%
%Khuloud Jaqaman, June 2012

%% Input

if nargin < 2
    diffModeDividerStruct = [];
end

if ~isempty(diffModeDividerStruct)
    
    %get trajectory lengths and positional standard deviations with diffusion
    %mode dividers
    trajLengthWithDivider = diffModeDividerStruct.trajLength;
    coordStdWithDivider = diffModeDividerStruct.coordStd;
    
    %thus get minimum trajectory length that can be classified
    minLength = min(trajLengthWithDivider);
    
    %also get maximum trajectory length with its own divider
    maxLength = max(trajLengthWithDivider);
    
    %get number of diffusion mode dividers (= number of modes - 1)
    dividerValues = diffModeDividerStruct.divider;
    numModeDiv = size(dividerValues,3);
    
else
    
    %take minimum trajectory length for which a diffusion coefficient can
    %be calculated as 5
    minLength = 5;
    
end

%get number of tracks
numTracks = length(tracksFinal);

%% Diffusion mode classification

%reserve memory for output parameters
diffModeAnalysisRes = repmat(struct('diffMode',[],'diffCoef',[],'msdF2F',[],...
    'meanPosStd',[]),numTracks,1);

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get current track's coordinates
    trackCoordCurrent = tracksFinal(iTrack).tracksCoordAmpCG;
    xCoord = trackCoordCurrent(:,1:8:end);
    yCoord = trackCoordCurrent(:,2:8:end);
    
    %calculate current track's displacements along x and y
    xCoordDelta = diff(xCoord,[],2);
    yCoordDelta = diff(yCoord,[],2);
    
    %calculate effective trajectory length from the number of available
    %displacements (this takes care of gaps)
    numSeg = size(xCoordDelta,1);
    segLft = NaN(numSeg,1);
    for iSeg = 1 : size(xCoordDelta,1)
        segLft(iSeg) = length(find(~isnan(xCoordDelta(iSeg,:)))) + 1;
    end
    
    %determine which segments are of sufficient length to be classified
    indxGood = find(segLft >= minLength);
    indxBad  = setdiff(1:numSeg,indxGood);
    
    %calculate the mean square frame-to-frame displacement
    msdF2F = nanmean(xCoordDelta.^2+yCoordDelta.^2,2);
        
    %get current track's positional standard deviations
    xCoordStd = trackCoordCurrent(:,5:8:end);
    yCoordStd = trackCoordCurrent(:,6:8:end);
    
    %calculate mean positional variance per track
    %also determine which coordStd is relevant
    meanPosVar = nanmean([xCoordStd yCoordStd].^2,2);
    meanPosStd = sqrt(meanPosVar);

    %calculate diffusion coefficient
    diffCoefCurrent = msdF2F/4 - meanPosVar;
    diffCoefCurrent(indxBad) = NaN;
    
    %store diffusion coefficient, mean square displacement and mean positional std in output variable
    diffModeAnalysisRes(iTrack).diffCoef = diffCoefCurrent;
    diffModeAnalysisRes(iTrack).msdF2F = msdF2F;
    diffModeAnalysisRes(iTrack).meanPosStd = meanPosStd;

    if ~isempty(diffModeDividerStruct)
        
        %determine diffusion mode
        diffModeCurrent = NaN(numSeg,1);
        for iSeg = indxGood'
            iLength = min(segLft(iSeg)-minLength+1,maxLength-minLength+1);
            tmp = abs(coordStdWithDivider-meanPosStd(iSeg));
            iStd = find(tmp==min(tmp));
            tmp = find(dividerValues(iLength,iStd,:)>diffCoefCurrent(iSeg),1,'first'); %#ok<FNDSB>
            if isempty(tmp)
                tmp = numModeDiv + 1;
            end
            diffModeCurrent(iSeg) = tmp;
        end
        
        %store diffusion mode in output variable
        diffModeAnalysisRes(iTrack).diffMode = diffModeCurrent;
        
    end
    
end

%call code to summarize diffusion analysis results
%store in .summary field of first track - rest stay empty
if isstruct(tracksFinal)
    probDim = 2;
    extractType = 1;
    [probMotionMode,modeMotionChar] = summarizeDiffModeRes(tracksFinal,...
        diffModeAnalysisRes,numModeDiv+1,minLength,probDim,extractType);
    diffModeAnalysisRes(1).summary.probMotionMode = probMotionMode;
    diffModeAnalysisRes(1).summary.modeMotionChar = modeMotionChar;
end

%% ~~~ the end ~~~

