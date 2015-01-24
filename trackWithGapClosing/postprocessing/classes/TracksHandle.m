classdef TracksHandle < Tracks
% TracksHandle is a Tracks implementation optimized for serving the new
% properties such as X, Y, Z while also providing backwards-compatability
% with the tracksFinal struct
%
% See also Tracks
%
% Mark Kittisopikul, January 2015
    properties
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
        % numSegments x numFrames matrix of indices. See class description.
        tracksFeatIndxCG
    end
    properties (Dependent = true)
        % numSegments x 8 x numFrames matrix of coordinates and amplitudes. See class description.
        tracksCoordAmpCG3D
        % 2D matrix, corresponds to tracksCoordAmpCG3D(:,:)
        tracksCoordAmpCG
        % numEvents x 4 matrix. See class description.
        seqOfEvents
        % Number of segments in each compound track
        % see also getNumSegments
        numSegments
        % Number of frames in which each compound track exists
        numFrames
    end
    properties (Access = protected, Transient)
        cache
    end
    methods
        function obj = TracksHandle(s)
            % Takes a tracksFinal structure from trackCloseGapsKalman
            if(nargin ~= 0)
                if(~isstruct(s))
                    s = convertMat2Struct2(s);
                end
                obj(numel(s)) = TracksHandle();
                [obj.tracksFeatIndxCG] = deal(s.tracksFeatIndxCG);
                [obj.tracksCoordAmpCG] = deal(s.tracksCoordAmpCG);
                [obj.seqOfEvents] = deal(s.seqOfEvents);
            end
        end
        function set.tracksFeatIndxCG(obj,tracksFeatIndxCG)
            obj.tracksFeatIndxCG = tracksFeatIndxCG;
        end
        function set.tracksCoordAmpCG(obj,tracksCoordAmpCG)
            threeD = reshape(tracksCoordAmpCG,size(tracksCoordAmpCG,1),8,[]);
            obj.tracksCoordAmpCG3D = threeD;
        end
        function tracksCoordAmpCG = get.tracksCoordAmpCG(obj)
            tracksCoordAmpCG = obj.tracksCoordAmpCG3D(:,:);
        end
        function set.tracksCoordAmpCG3D(obj,tracksCoordAmpCG3D)
            obj.X  = tracksCoordAmpCG3D(:,1,:);
            obj.Y  = tracksCoordAmpCG3D(:,2,:);
            obj.Z  = tracksCoordAmpCG3D(:,3,:);
            obj.A  = tracksCoordAmpCG3D(:,4,:);
            obj.dX = tracksCoordAmpCG3D(:,5,:);
            obj.dY = tracksCoordAmpCG3D(:,6,:);
            obj.dZ = tracksCoordAmpCG3D(:,7,:);
            obj.dA = tracksCoordAmpCG3D(:,8,:);
            
            obj.X  = obj.X(:,:);
            obj.Y  = obj.Y(:,:);
            obj.Z  = obj.Z(:,:);
            obj.A  = obj.A(:,:);
            obj.dX = obj.dX(:,:);
            obj.dY = obj.dY(:,:);
            obj.dZ = obj.dZ(:,:);
            obj.dA = obj.dA(:,:);
        end
        function tracksCoordAmpCG3D = get.tracksCoordAmpCG3D(obj)
            tracksCoordAmpCG3D = zeros(obj.numSegments,8,obj.numFrames);
            tracksCoordAmpCG3D(:,1,:) = obj.X;
            tracksCoordAmpCG3D(:,2,:) = obj.Y;
            tracksCoordAmpCG3D(:,3,:) = obj.Z;
            tracksCoordAmpCG3D(:,4,:) = obj.A;
            tracksCoordAmpCG3D(:,5,:) = obj.dX;
            tracksCoordAmpCG3D(:,6,:) = obj.dY;
            tracksCoordAmpCG3D(:,7,:) = obj.dZ;
            tracksCoordAmpCG3D(:,8,:) = obj.dA;
        end
        function set.seqOfEvents(obj,seqOfEvents)
            init = zeros(obj.numSegments,1);
            
            obj.segmentStartFrame = init;
            startIdx = seqOfEvents(:,2) == 1;
            obj.segmentStartFrame(seqOfEvents(startIdx,3)) = seqOfEvents(startIdx,1);
            
            obj.segmentEndFrame = init;
            endIdx = seqOfEvents(:,2) == 2;
            obj.segmentEndFrame(seqOfEvents(endIdx,3)) = seqOfEvents(endIdx,1);
            endIdx = endIdx & ~isnan(seqOfEvents(:,4));
            obj.segmentEndFrame(seqOfEvents(endIdx,3)) = seqOfEvents(endIdx,1) - 1;
            
            obj.parentSegment = init;
            idx = seqOfEvents(:,2) == 1;
            obj.parentSegment(seqOfEvents(idx,3)) = seqOfEvents(idx,4);
            
            obj.spouseSegment = init;
            idx = seqOfEvents(:,2) == 2;
            obj.spouseSegment(seqOfEvents(idx,3)) = seqOfEvents(idx,4);
            
            obj.startFrame = seqOfEvents(1,1);
            obj.endFrame = seqOfEvents(end,1);
        end
        function seqOfEvents = get.seqOfEvents(obj)
            seqOfEvents = zeros(obj.numSegments*2,4);
            
            startIdx = 1:obj.numSegments;
            seqOfEvents(startIdx,1) = obj.segmentStartFrame;
            seqOfEvents(startIdx,2) = 1;
            seqOfEvents(startIdx,3) = 1:obj.numSegments;
            seqOfEvents(startIdx,4) = obj.parentSegment;
            
            % For segments that merge, their true endFrames are offset by 1
            E = obj.segmentEndFrame;
            segmentMerged = ~isnan(obj.spouseSegment);
            E(segmentMerged) = E(segmentMerged) + 1;
            
            endIdx = (1:obj.numSegments) + obj.numSegments;
            seqOfEvents(endIdx,1) = E;
            seqOfEvents(endIdx,2) = 2;
            seqOfEvents(endIdx,3) = 1:obj.numSegments;
            seqOfEvents(endIdx,4) = obj.spouseSegment;
            
            seqOfEvents = sortrows(seqOfEvents);
        end
        function resetCache(obj)
            c = struct();
            c.tracksCoordAmpCG = [];
            obj.cache = c;
        end
        function N = get.numSegments(obj)
            N = size(obj.X,1);
        end
        function N = get.numFrames(obj)
            N = size(obj.X,2);
        end
    end
end