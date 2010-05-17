function [behaviorInWindows,tracksInWindows] = particleBehaviorInWindows(tracksFinal,...
    winPositions,winFrames,diffAnalysisRes,minLength)
%PARTICLEBEHAVIORINWINDOWS averages single particle behavior in windows based on cell edge segmentation
%
%SYNOPSIS 
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       winPositions   : The window edges for all time points, as output by
%                        Hunter's old windowing function.
%       winFrames      : The frames at which there are windows.
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%                        Optional. If not input but needed, it will be
%                        calculated within the code.
%       minLength      : Minimum length of a trajectory to include in
%                        analysis.
%                        Optional. Default: 5.
%
%OUTPUT behaviorInWindows: 4-D matrix storing average particle behavior
%                        property in each window. The dimensions are:
%                        1 - window index perpendicular to cell edge.
%                        2 - window index parallalel to cell edge.
%                        3 - window index in time.
%                        4 - property being averaged. They are the
%                        following:
%                        From track diffusion analysis ...
%                           (1) Fraction of unclassified tracks.
%                           (2) Among the classified tracks, fraction of confined.
%                           (3) //, fraction of Brownian.
%                           (4) //, fraction of directed.
%                           (5) Average diffusion coefficient.
%                           (6) Average confinement radius.
%                       From the tracks directly ...
%                           (7) Average frame-to-frame displacement.
%                           (8) Single particle density.
%                       Not implemented yet ...
%                           (9) Merge/split events.
%
%REMARKS This code is designed for experiments where the particle
%        trajectories are sampled much more frequently than the cell edge.
%        It assumes that particle lifetimes are much shorter than the time
%        between cell edge frames.
%
%        For a different scenario where particle lifetimes are longer than
%        the time between cell edge frames, the tracks should not be
%        grouped like this. Rather, each track should get divided into
%        several segments corresponding to the times between cell edge
%        frames and each track segment should be analyzed separately.
%        Something like that.
%
%Khuloud Jaqaman, May 2010

%% Input

if nargin < 3
    disp('--particleBehaviorInWindows: Incorrect number of input arguments!');
    return
end

if nargin < 4 || isempty(diffAnalysisRes)
    diffAnalysisRes = trackDiffusionAnalysis1(tracksFinal,1,2,0,0.05);
end

if nargin < 5 || isempty(minLength)
    minLength = 5;
end

%get the number of windows in each dimension
[numWinPerp,numWinPara,numWinFrames] = size(winPositions);

%define number of properties to be calculated
numProperties = 8;

%determine the number of SPT frames between window frames
numSPTFrames = mean(diff(winFrames));

%% Trajectory pre-processing

%keep only trajectories longer than minLength
criteria.lifeTime.min = minLength;
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx,:);
diffAnalysisRes = diffAnalysisRes(indx);

%convert tracksFinal into matrix if it's a structure
inputStructure = tracksFinal;
if isstruct(tracksFinal)
    clear tracksFinal
    tracksFinal = convStruct2MatIgnoreMS(inputStructure);
end

%% Particle behavior pre-processing

%From diffusion analysis ...

%get trajectory classifications
trajClass = vertcat(diffAnalysisRes.classification);
trajClass = trajClass(:,2);

%get diffusion coefficients
diffCoefGen = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');

%get confinement radii
confRad = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius(:,1)');

%From tracks directly ...

%extract the x- and y-coordinates from the track matrix
xCoord = tracksFinal(:,1:8:end);
yCoord = tracksFinal(:,2:8:end);

%calculate the average frame-to-frame displacemt per track
frame2frameDisp = nanmean( sqrt( diff(xCoord,[],2).^2 + diff(yCoord,[],2).^2 ) ,2);

%get the lifetime of each track
trackLft = getTrackSEL(tracksFinal);
trackLft = trackLft(:,3);

%% Window pre-processing

%initialize array storing window sizes
winSize = NaN(numWinPerp,numWinPara,numWinFrames-1);

%go over all windows and get their sizes
for iFrame = 1 : numWinFrames-1
    for iPara = 1 : numWinPara
        for iPerp = 1 : numWinPerp
            
            %get the window boundaries
            winX = [winPositions(iPerp,iPara,iFrame).outerBorder(1,:) ...
                winPositions(iPerp,iPara,iFrame).innerBorder(1,end:-1:1)]';
            winY = [winPositions(iPerp,iPara,iFrame).outerBorder(2,:) ...
                winPositions(iPerp,iPara,iFrame).innerBorder(2,end:-1:1)]';
            
            %calculate the window size
            winSize(iPerp,iPara,iFrame) = polyarea(winX,winY);
            
        end
    end
end

%% Calculate property values per window

%divide the trajectories among the windows
tracksInWindows = assignTracks2Windows(tracksFinal,winPositions,winFrames);

%initialize output variables
behaviorInWindows = NaN(numWinPerp,numWinPara,numWinFrames-1,numProperties);

%go over all windows and calculate particle properties in each
for iFrame = 1 : numWinFrames-1
    for iPara = 1 : numWinPara
        for iPerp = 1 : numWinPerp
            
            %get the tracks belonging to this window
            tracksCurrent = tracksInWindows{iPerp,iPara,iFrame};
            numTracksCurrent = length(tracksCurrent);
            
            %if there are tracks in this window ...
            if numTracksCurrent ~= 0
                
                %calculate the fraction of tracks in each motion category
                trajClassCurrent = trajClass(tracksCurrent);
                fracUnClass = length(find(isnan(trajClassCurrent)))/numTracksCurrent; %unclassified tracks
                trajClassCurrent = trajClassCurrent(~isnan(trajClassCurrent)); %classified tracks
                if ~isempty(trajClassCurrent)
                    fracClass = hist(trajClassCurrent,(1:3))/length(trajClassCurrent);
                else
                    fracClass = NaN(1,3);
                end
                
                %calculate the average diffusion coefficient
                diffCoefAve = nanmean(diffCoefGen(tracksCurrent));
                
                %calculate the average confinement radius
                confRadAve = nanmean(confRad(tracksCurrent));
                
                %calculate the average frame-to-frame displacement
                f2fDispAve = nanmean(frame2frameDisp(tracksCurrent));
                
                %calculate the average particle density
                partDensity = sum(trackLft(tracksCurrent)) / ...
                    winSize(iPerp,iPara,iFrame) / numSPTFrames;
                
                %store all properties in output cell array
                behaviorInWindows(iPerp,iPara,iFrame,:) = [fracUnClass ...
                    fracClass diffCoefAve confRadAve f2fDispAve partDensity];
            end
            
        end
    end
end

