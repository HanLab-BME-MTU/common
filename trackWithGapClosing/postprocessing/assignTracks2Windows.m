function tracksInWindows = assignTracks2Windows(tracksFinal,winPositions,...
    winFrames)
%ASSIGNTRACKS2WINDOWS groups tracks into spatial and temporal windows derived from the cell edge
%
%SYNOPSIS tracksInWindows = assignTracks2Windows(tracksFinal,winPositions,...
%    winFrames)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       winPositions   : The window edges for all time points, as output by
%                        Hunter's old windowing function.
%       winFrames      : The frames at which there are windows.
%
%OUTPUT tracksInWindows: Cell array of dimensions (number of bands) x
%                        (number of windows) x (number of window frames -1)
%                        storing the track indices that fall in each window
%                        in each frame.
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
    disp('--assignTracks2Windows: Incorrect number of input arguments!');
    return
end

%get number of tracks
numTracks = size(tracksFinal,1);

%get number of windows along the edge and perpendicular to it, and number
%of frames that have windows
[numWinPerp,numWinPara,numWinFrames] = size(winPositions);

%% Pre-processing

%get track start and end times
trackSEL = getTrackSEL(tracksFinal);

%calculate the "average" time at which a track exists
%this will be used to assign tracks to time windows
trackTimeMean = mean(trackSEL(:,1:2),2);

%get average track positions
if isstruct(tracksFinal) %if compound tracks
    
    [xCoordMean,yCoordMean] = deal(NaN(numTracks,1));
    for iTrack = 1 : numTracks
        xCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
        yCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
        xCoordMean(iTrack) = nanmean(xCoordAll(:));
        yCoordMean(iTrack) = nanmean(yCoordAll(:));
    end
    
else %if simple tracks
    
    %extract x- and y-coordinates
    xCoordAll = tracksFinal(:,1:8:end);
    yCoordAll = tracksFinal(:,2:8:end);
    
    %calculate the average coordinates
    xCoordMean = nanmean(xCoordAll,2);
    yCoordMean = nanmean(yCoordAll,2);
    
end

%% Track assignment into windows

%initialize cell array storing the grouping of tracks into windows for
%each frame range
tracksInWindows = cell(numWinPerp,numWinPara,numWinFrames-1);

%go over all window frames
for iWinFrame = 1 : numWinFrames - 1
    
    %get current frame number and next frame number
    minFrame = winFrames(iWinFrame);
    maxFrame = winFrames(iWinFrame + 1);
    
    %find tracks whose "average" time is in this frame range
    indxFrameRange = find(trackTimeMean>=minFrame & trackTimeMean<maxFrame);
    
    %get the "average" positions of these tracks
    xCoordMeanFR = xCoordMean(indxFrameRange);
    yCoordMeanFR = yCoordMean(indxFrameRange);
    
    %go over the windows in this frame
    for iPara = 1 : numWinPara
        for iPerp = 1 : numWinPerp
            
            %get the window boundaries
            winX = [winPositions(iPerp,iPara,iWinFrame).outerBorder(1,:) ...
                winPositions(iPerp,iPara,iWinFrame).innerBorder(1,end:-1:1)]';
            winY = [winPositions(iPerp,iPara,iWinFrame).outerBorder(2,:) ...
                winPositions(iPerp,iPara,iWinFrame).innerBorder(2,end:-1:1)]';

            %find the tracks whose "average" position lies in this window
            indxWin = inpolygon(xCoordMeanFR,yCoordMeanFR,winX,winY);
            
            %map back to original track indices
            indxWin = indxFrameRange(indxWin);
            
            %store track indices in cell array
            tracksInWindows{iPerp,iPara,iWinFrame} = indxWin;
            
        end
    end
    
end
