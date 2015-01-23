classdef Tracks < handle &  matlab.mixin.Copyable
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
    properties
        % numSegments x numFrames matrix of indices. See class description.
        tracksFeatIndxCG
        % numSegments x 8 x numFrames matrix of coordinates and amplitudes. See class description.
        tracksCoordAmpCG3D
        % numEvents x 4 matrix. See class description.
        seqOfEvents
    end
    properties (Dependent = true)
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
    end
    properties (Access = protected, Transient)
        cache
    end
    methods
        function obj = Tracks(s)
            % Takes a tracksFinal structure from trackCloseGapsKalman
            if(nargin ~= 0)
                if(~isstruct(s))
                    s = convertMat2Struct2(s);
                end
                obj(numel(s)) = Tracks();
                [obj.tracksFeatIndxCG] = deal(s.tracksFeatIndxCG);

                threeD = cellfun(@(c) reshape(c,size(c,1),8,[]), ...
                    {s.tracksCoordAmpCG},'UniformOutput',false);
                
                [obj.tracksCoordAmpCG3D] = deal(threeD{:});
                [obj.seqOfEvents] = deal(s.seqOfEvents);
            end
        end
        function resetCache(obj)
            c = struct();
            obj.cache = c;
        end
        function set.tracksFeatIndxCG(obj,tracksFeatIndxCG)
            obj.tracksFeatIndxCG = tracksFeatIndxCG;
            obj.resetCache();
        end
        function set.tracksCoordAmpCG3D(obj,tracksCoordAmpCG3D)
            obj.tracksCoordAmpCG3D = tracksCoordAmpCG3D;
            obj.resetCache();
        end
        function set.seqOfEvents(obj,seqOfEvents)
            obj.seqOfEvents = seqOfEvents;
            obj.resetCache();
        end
        function tracksCoordAmpCG = get.tracksCoordAmpCG(obj)
%             tracksCoordAmpCG = obj.tracksCoordAmpCG3D(:,:);
            tracksCoordAmpCG = reshape(obj.tracksCoordAmpCG3D,size(obj.tracksCoordAmpCG3D,1),[]);;
        end
        function set.tracksCoordAmpCG(obj,tracksCoordAmpCG)
            obj.tracksCoordAmpCG3D = ... 
                reshape(tracksCoordAmpCG,size(tracksCoordAmpCG,1),8,[]);
        end
        function X = get.X(obj)
            X = obj.tracksCoordAmpCG3D(:,1,:);
            X = X(:,:);
        end
        function Y = get.Y(obj)
            Y = obj.tracksCoordAmpCG3D(:,2,:);
            Y = Y(:,:);
        end
        function Z = get.Z(obj)
            Z = obj.tracksCoordAmpCG3D(:,3,:);
            Z = Z(:,:);
        end
        function A = get.A(obj)
            A = obj.tracksCoordAmpCG3D(:,4,:);
            A = A(:,:);
        end
        function dX = get.dX(obj)
            dX = obj.tracksCoordAmpCG3D(:,5,:);
            dX = dX(:,:);
        end
        function dY = get.dY(obj)
            dY = obj.tracksCoordAmpCG3D(:,6,:);
            dY = dY(:,:);
        end
        function dZ = get.dZ(obj)
            dZ = obj.tracksCoordAmpCG3D(:,7,:);
            dZ = dZ(:,:);
        end
        function dA = get.dA(obj)
            dA = obj.tracksCoordAmpCG3D(:,8,:);
            dA = dA(:,:);
        end
        function S = get.segmentStartFrame(obj)
            S = zeros(obj.numSegments,1);
            startIdx = obj.seqOfEvents(:,2) == 1;
            S(obj.seqOfEvents(startIdx,3)) = obj.seqOfEvents(startIdx,1);
        end
        function E = get.segmentEndFrame(obj)
            E = zeros(obj.numSegments,1);
            endIdx = obj.seqOfEvents(:,2) == 2;
            E(obj.seqOfEvents(endIdx,3)) = obj.seqOfEvents(endIdx,1);
            endIdx = endIdx & ~isnan(obj.seqOfEvents(:,4));
            E(obj.seqOfEvents(endIdx,3)) = obj.seqOfEvents(endIdx,1) - 1;
        end
        function O = get.parentSegment(obj)
            O = zeros(obj.numSegments,1);
            idx = obj.seqOfEvents(:,2) == 1;
            O(obj.seqOfEvents(idx,3)) = obj.seqOfEvents(idx,4);
        end
        function D = get.spouseSegment(obj)
            D = zeros(obj.numSegments,1);
            idx = obj.seqOfEvents(:,2) == 2;
            D(obj.seqOfEvents(idx,3)) = obj.seqOfEvents(idx,4);
        end
        function startTime = get.startFrame(obj)
            startTime = obj.seqOfEvents(1,1);
        end
        function endTime = get.endFrame(obj)
            endTime = obj.seqOfEvents(end,1);
        end
        out = textGraph(obj)
        disp(obj)
        plot(obj,varargin)
        m = getMatrix(obj,varargin)
        s = getSparse(obj,varargin)
        n = numTimePoints(obj);
        n = numSegments(obj,idx);
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
