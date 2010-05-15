function particleBehaviorVsEdgeMotion(tracksFinal,winPositions,winFrames,...
    protrusionWin,diffAnalysisRes,minLength)
%PARTICLEBEHAVIORVSEDGEMOTION looks for correlation between particle behavior and cell edge protrusion activity
%
%SYNOPSIS 
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       winPositions   : The window edges for all time points, as output by
%                        Hunter's old windowing function.
%       winFrames      : The frames at which there are windows.
%       protrusionWin  : Average protrusion vector per window (from
%                        Hunter).
%       diffAnalysisRes: Output of trackDiffusionAnalysis1.
%                        Optional. If not input but needed, it will be
%                        calculated within the code.
%       minLength      : Minimum length of a trajectory to include in
%                        analysis.
%                        Optional. Default: 20.
%
%OUTPUT 
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

if nargin < 4
    disp('--particleBehaviorVsEdgeMotion: Incorrect number of input arguments!');
    return
end

if nargin < 5 || isempty(diffAnalysisRes)
    diffAnalysisRes = [];
end

if nargin < 6 || isempty(minLength)
    minLength = 20;
end


