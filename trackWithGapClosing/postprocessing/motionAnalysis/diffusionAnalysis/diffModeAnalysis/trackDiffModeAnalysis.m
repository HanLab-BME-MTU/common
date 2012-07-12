function diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDivider,minLength)
%TRACKDIFFMODEANALYSIS classifies tracks into diffusion modes
%
%SYNOPSIS diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDivider,minLength)
%
%INPUT  tracksFinal : Output of trackCloseGapsKalman.
%       diffModeDivider: Vector of values dividing between the different
%                     diffusion modes. If tehre are N modes, this will be a
%                     vector of N-1 entries.
%       minLength   : Minimum length of a track to be included in analysis.
%                     Optional. Default: 5 frames.
%
%OUTPUT diffModeAnalysisRes : Structure array with the following fields per
%                         track:
%           .diffMode: Diffusion mode.
%           .diffCoef: Diffusion coefficient as calculated from the mean
%           square frame-to-frame displacement.
%
%REMARKS Code is written for 2D case only, but can be generalized to 3D.
%
%Khuloud Jaqaman, June 2012

%% Input

if nargin < 3 || isempty(minLength)
    minLength = 5;
end

%get number of tracks
numTracks = length(tracksFinal);

%get number of diffusion mode divider (= number of modes - 1)
numModeDiv = length(diffModeDivider);

%% Diffusion mode classification

%reserve memory for output parameters
diffModeAnalysisRes = repmat(struct('diffMode',[],'diffCoef',[]),numTracks,1);

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get current track's coordinates
    trackCoordCurrent = tracksFinal(iTrack).tracksCoordAmpCG;
    xCoord = trackCoordCurrent(:,1:8:end);
    yCoord = trackCoordCurrent(:,2:8:end);
    
    %calculate current track's displacements along x and y
    xCoordDelta = diff(xCoord,[],2);
    yCoordDelta = diff(yCoord,[],2);
    
    %calculate the mean square frame-to-frame displacement
    msdF2F = nanmean(xCoordDelta.^2+yCoordDelta.^2,2);
        
    %get current track's positional standard deviations
    xCoordStd = trackCoordCurrent(:,5:8:end);
    yCoordStd = trackCoordCurrent(:,6:8:end);
    
    %calculate mean positional variance per track
    meanPosVar = nanmean([xCoordStd yCoordStd].^2,2);

    %determine which segments are of sufficient length
    segLft = getTrackSEL(trackCoordCurrent);
    numSeg = size(segLft,1);
    indxGood = find(segLft(:,3) >= minLength);
    indxBad  = setdiff(1:numSeg,indxGood);
    
    %calculate diffusion coefficient
    diffCoefCurrent = msdF2F/4 - meanPosVar;
    diffCoefCurrent(indxBad) = NaN;
    
    %determine diffusion mode
    diffModeCurrent = diffCoefCurrent;
    diffModeCurrent(diffCoefCurrent<diffModeDivider(1)) = 1;
    for iDiv = 1 : numModeDiv
        diffModeCurrent(diffCoefCurrent>diffModeDivider(iDiv)) = iDiv + 1;
    end
    
    %store results in output variable
    diffModeAnalysisRes(iTrack).diffCoef = diffCoefCurrent;
    diffModeAnalysisRes(iTrack).diffMode = diffModeCurrent;
    
end


%% ~~~ the end ~~~

