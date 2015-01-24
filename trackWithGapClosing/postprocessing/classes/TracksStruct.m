classdef TracksStruct < Tracks
% TracksStruct is a Tracks implementation that is a close extension of the
% original struct implementation for backwards compatibility while also
% providing the properties and methods of the Tracks interface class.
%
% Emphasis: Backwards-compatability, fast-instantiation
%
% See also Tracks
%
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
        % Number of segments in each compound track
        % see also getNumSegments
        numSegments
        % Number of frames in which each compound track exists
        numFrames
    end
    methods
        function obj = TracksStruct(s)
            % Takes a tracksFinal structure from trackCloseGapsKalman
            if(nargin ~= 0)
                if(~isstruct(s))
                    s = convertMat2Struct2(s);
                end
                obj(numel(s)) = TracksStruct();
                [obj.tracksFeatIndxCG] = deal(s.tracksFeatIndxCG);
                [obj.tracksCoordAmpCG] = deal(s.tracksCoordAmpCG);
                [obj.seqOfEvents] = deal(s.seqOfEvents);
            end
        end
        function tracksCoordAmpCG = get.tracksCoordAmpCG(obj)
            if(isempty(obj.cache.tracksCoordAmpCG))
            tracksCoordAmpCG = obj.tracksCoordAmpCG3D(:,:);
%                 tracksCoordAmpCG = reshape(obj.tracksCoordAmpCG3D,size(obj.tracksCoordAmpCG3D,1),[]);
                obj.cache.tracksCoordAmpCG = tracksCoordAmpCG;
            else
                tracksCoordAmpCG = obj.cache.tracksCoordAmpCG;
            end
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
        function N = get.numSegments(obj)
            N = size(obj.tracksCoordAmpCG3D,1);
        end
        function N = get.numFrames(obj)
            N = size(obj.tracksCoordAmpCG3D,3);
        end
    end
end