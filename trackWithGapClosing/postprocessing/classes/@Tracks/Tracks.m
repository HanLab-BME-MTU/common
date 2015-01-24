classdef (Abstract = true) Tracks < handle  & matlab.mixin.Copyable
% Tracks is a class that encapsulates the tracksFinal output of
% TrackingProcess. Each individual Tracks object consists of multiple
% segments that may merge and split. An array of Tracks objects represents
% individual tracks in the same dataset.
%
% This class extends the struct that is the output of trackCloseGapsKalman.
%
% See also trackObj , plotTracks2D, 
% common/trackWithGapClosing/postprocessing/misc
% common/trackWithGapClosing/postprocessing/motionAnalysis
% common/trackWithGapClosing/postprocessing/HuetImplementation
%
% Properties
%
% .tracksFeatIndxCG  : Indices of the detections in a movieInfo
%                      structure
%
% .tracksCoordAmpCG3D: The positions and amplitudes of the tracked
%                      features, after gap closing. 
%                               
%                      This 3D matrix has dimensions as follows:
%                      numSegments x 8 x numFrames
%                      * numSegments (rows) = number of track segments in compound
%                        track.
%                      * 8 corresponds to [x, y, z, a, dx, dy, dz, da]
%                      * numFrames = the number of time points in each
%                        segment. Column number is relative to
%                        .seqOfEvents(1) or .startFrame
%
%                      NaN indicates frames where track segments do not exist.
%
%                      Due to linear indexing this is mostly backwards
%                      compatible with the tracksCoordAmpCG property below
%                      and in the original structure when indexed.
%                      e.g. tracksCoordAmpCG3D(:,1:8:end) still accesses the
%                      x coordinate
%                      The incompatibility is when this is referenced
%                      without indexing.
%
% .seqOfEvents      : Matrix with number of rows equal to number
%                     of events happening in a track and 4
%                     columns:
%                     1st: Frame where event happens;
%                     2nd: 1 - start of track, 2 - end of track;
%                     3rd: Index of track segment that ends or starts;
%                     4th: NaN - start is a birth and end is a death,
%                          number - start is due to a split,
%                                   end is due to a merge
%                                   number is the index of track segment
%                                   for the merge/split.
% Dependent Properties
%
% X, Y, Z, A, dX, dY, dZ, dA : xyz coordinates, amplitude and uncertanties
% segmentStartFrame          : Absolute frame where each segment starts
% segmentEndFrame            : Absolute frame where each segment ends
% parentSegment              : Segment from which each segment split.
%                              NaN if the segment originated independently
% spouseSegment              : Segment into which the segment merged
%                              NaN if the segment never merged

% Mark Kittisopikul, January 2015
    properties (Abstract = true)
        % numSegments x numFrames matrix of indices. See class description.
        tracksFeatIndxCG
        % numSegments x 8 x numFrames matrix of coordinates and amplitudes. See class description.
        tracksCoordAmpCG3D
        % numEvents x 4 matrix. See class description.
        seqOfEvents
        % 2D matrix, corresponds to tracksCoordAmpCG3D(:,:)
        tracksCoordAmpCG
        % X coordinates as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,1:8:end)
        % Column is relative to the startFrame
        X
        % Y coordinates as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,2:8:end)
        % Column is relative to the startFrame
        Y
        % Z coordinates as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,3:8:end)
        % Column is relative to the startFrame
        Z
        % Amplitude as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,4:8:end)
        % Column is relative to the startFrame
        A
        % X uncertainty as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,5:8:end)
        % Column is relative to the startFrame
        dX
        % Y uncertainty as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,6:8:end)
        % Column is relative to the startFrame
        dY
        % Z uncertainty as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,7:8:end)
        % Column is relative to the startFrame
        dZ
        % A uncertainty as nSeg x nFrame matrix = tracksCoordAmpCG3D(:,8:8:end)
        % Column is relative to the startFrame
        dA
        % Column vector for the absolute time point each segment started
        segmentStartFrame
        % Column vector for the absolute time point each segmented ended
        segmentEndFrame
        % Column vector for the segment from which the segment originated
        parentSegment
        % Column vector for the segment into which the segment merged
        spouseSegment
        % Scalar absolute frame at which the compound track starts
        startFrame
        % Scalar absolute at which the compound track ends
        endFrame
        % Number of segments in each compound track
        % see also getNumSegments
        numSegments
        % Number of frames in which each compound track exists
        numFrames
    end
    methods
        out = textGraph(obj)
        disp(obj)
        plot(obj,varargin)
        m = getMatrix(obj,varargin)
        s = getSparse(obj,varargin)
        n = numTimePoints(obj);
        t = totalSegments(obj)
        b = isstruct(obj)
        s = getStruct(obj)
        seqM = getSeqOfEventsMatrix(obj)
        msM = getMergeSplitMatrix(obj)
        idx = getMergeIdx(obj,msM)
        idx = getSplitIdx(obj,msM)
        [merge,split] = getMergeSplitXY(obj,matrix);
    end
end
