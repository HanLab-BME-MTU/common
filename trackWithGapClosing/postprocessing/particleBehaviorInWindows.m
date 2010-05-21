function [sptPropInWindow,tracksInWindow] = particleBehaviorInWindows(tracksFinal,...
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
%OUTPUT sptPropInWindow: Structure storing various single particle
%                        properties, with the following fields:
%               From the tracks directly ...
%           .spDensity      : Single particle density.
%           .f2fDisp        : Average frame-to-frame displacement.
%               From track diffusion analysis ...
%           .fracUnclass    : Fraction of unclassified tracks.
%           .fracConf       : Among the classified tracks, fraction of confined.
%           .fracBrown      : //, fraction of Brownian.
%           .fracDir        : //, fraction of directed.
%           .diffCoef       : Average diffusion coefficient.
%           .confRad        : Average confinement radius.
%                        Each field contains 2 sub-fields, "values" and
%                        "numPoints" storing the values of the properties
%                        and the number of points contributing to each
%                        value, respectively.
%                        Each of these fields is a 3-D matrix, of
%                        dimensions (number of bands) x (number of windows)
%                        x (number of window frames-1).
%       tracksInWindow : Output of assignTracks2Windows. Cell array of
%                        dimensions (number of bands) x (number of windows)
%                        x (number of window frames-1) storing the track 
%                        indices that fall in each window in each frame.
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

%determine the number of SPT frames between window frames
numSPTFrames = winFrames(2) - winFrames(1);

%% Trajectory pre-processing

%keep only trajectories longer than minLength
criteria.lifeTime.min = minLength;
indx = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx,:);
diffAnalysisRes = diffAnalysisRes(indx);

%get number of trajectories
numTracks = length(indx);

%divide the trajectories among the windows
tracksInWindow = assignTracks2Windows(tracksFinal,winPositions,winFrames,1);

%% Particle behavior pre-processing

%From diffusion analysis ...

%get trajectory classifications
trajClass = vertcat(diffAnalysisRes.classification);
trajClass = trajClass(:,2);

%get diffusion coefficients
diffCoefGen = catStruct(1,'diffAnalysisRes.fullDim.genDiffCoef(:,3)');

%get confinement radii
confRadAll = catStruct(1,'diffAnalysisRes.confRadInfo.confRadius(:,1)');

%From tracks directly ...

%get the average frame-to-frame displacement
if isstruct(tracksFinal) %if tracks are in structre format
    
    %reserve memory for frame-to-frame displacement
    frame2frameDisp = NaN(length(trajClass),1);
    
    %initialize global segment index
    iSeg = 0;
    
    %go over all compound tracks
    for iTrack = 1 : numTracks
        
        %get current track's coordinates
        trackCoordCurrent = tracksFinal(iTrack).tracksCoordAmpCG;
        xCoord = trackCoordCurrent(:,1:8:end);
        yCoord = trackCoordCurrent(:,2:8:end);
        
        %calculate average frame-to-frame displacement
        f2fDispCurrent = nanmean( sqrt( diff(xCoord,[],2).^2 + diff(yCoord,[],2).^2 ) ,2);
        
        %get number of segments in this compound track
        numSeg = length(f2fDispCurrent);
        
        %store in big vector
        frame2frameDisp(iSeg+1:iSeg+numSeg) = f2fDispCurrent;
        
        %update global segment index
        iSeg = iSeg + numSeg;
        
    end
    
    
else %if tracks are in matrix format
    
    %extract the x- and y-coordinates from the track matrix
    xCoord = tracksFinal(:,1:8:end);
    yCoord = tracksFinal(:,2:8:end);
    
    %calculate the average frame-to-frame displacemt per track
    frame2frameDisp = nanmean( sqrt( diff(xCoord,[],2).^2 + diff(yCoord,[],2).^2 ) ,2);
    
end

%get the lifetime of each track
trackLft = getTrackSEL(tracksFinal,1);
trackLft = trackLft(:,3);

%% Window pre-processing

%initialize array storing window sizes
winSize = NaN(numWinPerp,numWinPara,numWinFrames-1);

%go over all windows and get their sizes
for iFrame = 1 : numWinFrames-1
    for iPara = 1 : numWinPara
        for iPerp = 1 : numWinPerp
            
            %if this window has a finite size
            if ~isempty(winPositions(iPerp,iPara,iFrame).outerBorder) ...
                    && ~isempty(winPositions(iPerp,iPara,iFrame).innerBorder)
                
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
end

%% Calculate property values per window

%initialize output variables
[spDensity,f2fDisp,fracUnclass,fracConf,fracBrown,fracDir,diffCoef,...
    confRad] = deal(struct('values',NaN(numWinPerp,numWinPara,numWinFrames-1),...
    'numPoints',zeros(numWinPerp,numWinPara,numWinFrames-1)));

%go over all windows and calculate particle properties in each
for iFrame = 1 : numWinFrames-1
    for iPara = 1 : numWinPara
        for iPerp = 1 : numWinPerp
            
            %get the tracks belonging to this window
            tracksCurrent = tracksInWindow{iPerp,iPara,iFrame};
            numTracksCurrent = length(tracksCurrent);
            
            %if there are tracks in this window ...
            if numTracksCurrent ~= 0
                
                %From the tracks directly ...
                
                %calculate the average particle density
                spDensity.values(iPerp,iPara,iFrame) = sum(trackLft(tracksCurrent)) / ...
                    winSize(iPerp,iPara,iFrame) / numSPTFrames;
                spDensity.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %calculate the average frame-to-frame displacement
                f2fDisp.values(iPerp,iPara,iFrame) = nanmean(frame2frameDisp(tracksCurrent));
                f2fDisp.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;
                
                %From the diffusion analysis ...

                %calculate the fraction of tracks in each motion category
                trajClassCurrent = trajClass(tracksCurrent);
                fracUnclass.values(iPerp,iPara,iFrame) = ...
                    length(find(isnan(trajClassCurrent)))/numTracksCurrent; %unclassified tracks
                fracUnclass.numPoints(iPerp,iPara,iFrame) = numTracksCurrent;

                trajClassCurrent = trajClassCurrent(~isnan(trajClassCurrent)); %classified tracks
                numTracksClass = length(trajClassCurrent);
                if ~isempty(trajClassCurrent)
                    fracClass = hist(trajClassCurrent,(1:3))/length(trajClassCurrent);
                    fracConf.values(iPerp,iPara,iFrame) = fracClass(1);
                    fracConf.numPoints(iPerp,iPara,iFrame) = numTracksClass;
                    fracBrown.values(iPerp,iPara,iFrame) = fracClass(2);
                    fracBrown.numPoints(iPerp,iPara,iFrame) = numTracksClass;
                    fracDir.values(iPerp,iPara,iFrame) = fracClass(3);
                    fracDir.numPoints(iPerp,iPara,iFrame) = numTracksClass;
                end
                
                %calculate the average diffusion coefficient
                tmpVec = diffCoefGen(tracksCurrent);
                tmpVec = tmpVec(~isnan(tmpVec));
                diffCoef.values(iPerp,iPara,iFrame) = mean(tmpVec);
                diffCoef.numPoints(iPerp,iPara,iFrame) = length(tmpVec);
                
                %calculate the average confinement radius
                tmpVec = confRadAll(tracksCurrent);
                tmpVec = tmpVec(~isnan(tmpVec));
                confRad.values(iPerp,iPara,iFrame) = mean(tmpVec);
                confRad.numPoints(iPerp,iPara,iFrame) = length(tmpVec);

            end %(if numTracksCurrent ~= 0)
            
        end %(for iPerp = 1 : numWinPerp)
    end %(for iPara = 1 : numWinPara)
end %(for iFrame = 1 : numWinFrames-1)

%store all properties in output structure
sptPropInWindow = struct('spDensity',spDensity,'f2fDisp',f2fDisp,...
    'fracUnclass',fracUnclass,'fracConf',fracConf,'fracBrown',fracBrown,...
    'fracDir',fracDir,'diffCoef',diffCoef,'confRad',confRad);



