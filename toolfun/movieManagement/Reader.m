classdef  Reader < handle
    % Reader is an abstract class for reading image data
    %
    % Known subclasses include
    %
    % BioFormatsReader
    % TiffSeriesReader
    % HCSReader
    %
    % Generally methods take the channel as an argument
    %   which is usually added by similar calls in Channel
    %
    % Arguments are in the order C, T, Z
    %
    % Subclasses should consider consider overloading
    % protected methods ending with an underscore, '_'
    % such as loadImage_ and loadStack_.
    %
    % The public loadImage and loadStack methods contain
    % input validation and default parameters which should
    % be changed minimally.
    %
    
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
    end
    
    methods ( Access = public )
        function I = loadImage(obj, c, t, varargin)
        % loadImage reads a single image plane as a 2D, YX Matrix
        %
        % loadImage(c,t,z) reads the YX plane at (c,t,z)
        %
        % loadImage(c,t) reads the YX plane at (c,t,1)
        %
        % Note: Override if you need to overload or change how the input is
        % checked. Otherwise override loadImage_.
        %
        % Example:
        %   reader = movieData.getReader();
        %   I = reader.loadImage(1,1);
        %   imshow(I,[]);
        %
        
        % Backwards compatability before 2014/08/18:
        % Previously this function was abstract and therefore should be
        % overridden by all subclasses written prior to 2014/08/18
         
            ip = inputParser;
            ip.addRequired('c', ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC()));
            % Parse c first for validation check
            %  since getSizeT and getSizeZ may depend on a valid c.
            ip.parse(c);
            ip.addRequired('t', ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT(c)));
            ip.addOptional('z', 1, ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeZ(c)));
            ip.parse(c, t, varargin{:});
            
            I = obj.loadImage_(obj, c , t , ip.Results.z);
        end
        
        function I = loadStack(obj, c, t, varargin)
        % loadStack reads a Z-stack and returns a YXZ Matrix
        %
        %   loadStack(c,t,Z) returns a Z-stack for channel c and
        %   time frame t  with Z planes indicated by the vector Z
        %
        %   loadStack(c,t) returns loadStack(c,t, 1: getSizeZ())
        %
        % output: a 3D matrix with dimensions YXZ
        %
        % Note: Override only if need you need to overload or change how
        % the input is checked. Otherwise override loadStack_.
        %
        % Example:
        %   reader = movieData.getReader();
        %   % The following are all equivalent
        %   zStackMatrix = reader.loadStack(1,1,1:reader.getSizeZ());
        %   zStackMatrix = reader.loadStack(1,1);
        %   zStackMatrix = reader.loadStack(1);
        %
       
            % Input check
            ip = inputParser;
            ip.addRequired('c', ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC() ) );
            % Parse c first for validation check
            %  since getSizeT and getSizeZ may depend on a valid c.
            ip.parse(c);
            ip.addRequired('t', ...
                @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT() ) );
            ip.addOptional('z', 1 : obj.getSizeZ(), ...
                @(x) all(ismember(x, 1 : obj.getSizeZ() ) ) );
            ip.parse(c, t, varargin{:});
            t = ip.Results.t;
            z = ip.Results.z;
            
            I = obj.loadStack_(c,t,z);
        end
    end
    methods( Access = protected )
        function I = loadStack_(obj, c, t, z)
        % loadStack_ is a validation-free version called by loadStack
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
        
            % Load first image to get class and dimensions.
            first = obj.loadImage_( c , t , z(1) );
            % "I" will be a YXZ matrix.
            I = zeros( [ size(first) length(z) ] , class(first));
            I(:,:,1) = first;
            
            for zi = 2:length(z)
                I(:,:,zi) = obj.loadImage_( c , t , z(zi) );
            end
                
        end
    %end
    %methods( Abstract, Access = protected )
        % I = loadImage_(obj, c, t, z )
        function I = loadImage_(obj, c, t, z )
            % loadImage_ is a validation-free implementation called by
            %   loadImage (or will be in the near future)
            %
            % Forward to validation version for now.
            % This WILL be an abstract method in the future.

            I = obj.loadImage(c,t,z);
        end
    end
    
end
