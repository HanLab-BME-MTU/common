classdef  Reader < handle
    % Concrete implementation of MovieObject for a single movie
    
    properties
        sizeX
        sizeY
        sizeZ
        sizeC
        sizeT
        bitDepth
    end
    
    methods(Abstract)
        getSizeX(obj)
        getSizeY(obj)
        getSizeZ(obj)
        getSizeC(obj)
        getSizeT(obj)
        getBitDepth(obj)
        getImageFileNames(obj)
        getChannelNames(obj)
        loadImage(obj)
    end
    
    methods
        % loadStack
        %
        % Provides a generic implementation of loadStack.
        % Guarantees all readers will have a loadStack method if they have
        % a loadImage method.
        %
        % Override if there are backend optimizations.
        % This is likely slower because of repeated input validation.
        %
        % parameters are c, t, z with only c required
        % t defaults to 1
        % z defaults to 1:sizeZ
        %
        % output a 3D matrix, YXZ
        function I = loadStack(obj, c, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('c', ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC() ) );
            % Parse c first for validation check
            %  since getSizeT and getSizeZ may depend on a valid c.
            ip.parse(c);
            ip.addOptional('t', 1, ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT(c) ) );
            ip.addOptional('z', 1 : obj.getSizeZ(c), ...
                @(x) all(ismember(x, 1 : obj.getSizeZ(c) ) ) );
            ip.parse(c, varargin{:});
            t = ip.Results.t;
            z = ip.Results.z;

            % Load first image to get class and dimensions.
            first = obj.loadImage( c , t(1) , z(1) );
            % "I" will be a YXZ matrix.
            I = zeros( [ size(first) length(z) ] , 'like', first);
            I(:,:,1) = first;
            
            for zi = 2:length(z)
                I(:,:,zi) = obj.loadImage( c , t(1) , z(zi) );
            end
                
        end
    end
    
end
