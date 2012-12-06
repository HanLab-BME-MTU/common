classdef  OMEROReader < Reader
    % Concrete implementation of MovieObject for a single movie
    properties
        imageID
    end
    
    properties (Transient = true)
        session
        image
        pixels
    end
    
    methods
        %% Constructor
        function obj = OMEROReader(imageID, varargin)
            
            obj.imageID = imageID;
            if nargin > 1, obj.setSession(varargin{1}); end
        end
        
        
        %% Dimensions functions
        function sizeX = getSizeX(obj, varargin)
            sizeX = obj.getPixels().getSizeX.getValue;
        end
        
        function sizeY = getSizeY(obj, varargin)
            sizeY = obj.getPixels().getSizeY.getValue;
        end
        
        function sizeC = getSizeC(obj, varargin)
            sizeC = obj.getPixels().getSizeC.getValue;
        end
        
        function sizeZ = getSizeZ(obj, varargin)
            sizeZ = obj.getPixels().getSizeZ.getValue;
        end
        
        function sizeT = getSizeT(obj, varargin)
            sizeT = obj.getPixels().getSizeT.getValue;
        end
        
        function session = getSession(obj)
            % Check session is not empty
            assert(~isempty(obj.session), 'No session created');
            session =  obj.session;
        end
        
        function setSession(obj, session)
            % Check input
            ip = inputParser;
            ip.addRequired('session', @(x) isa(x, 'omero.api.ServiceFactoryPrxHelper'));
            ip.parse(session);
            
            obj.session = session;
        end
        
        function fileNames = getImageFileNames(obj, iChan, varargin)
            % Generate image file names
            basename = sprintf('Image%g_c%d_t', obj.imageID, iChan);
            fileNames = arrayfun(@(t) [basename num2str(t, ['%0' num2str(floor(log10(obj.getSizeT))+1) '.f']) '.tif'],...
                1:obj.getSizeT,'Unif',false);
            
        end
        
        function chanNames = getChannelNames(obj, iChan)
            chanNames = arrayfun(@(x) ['Image ' num2str(obj.imageID) ...
                ': Channel ' num2str(x)], iChan, 'UniformOutput', false);
        end
        
        %% Image loading function
        function I = loadImage(obj, iChan, iFrame)
            % Test session integrity
            store = obj.getSession().createRawPixelsStore();
            store.setPixelsId(obj.getPixels().getId().getValue(), false);
            I = zeros([obj.getSizeY(), obj.getSizeX(), numel(iFrame)]);
            for i=1:numel(iFrame),
                plane = store.getPlane(0, iChan-1, iFrame(i)-1);
                I(:,:,i)=double(toMatrix(plane, obj.getPixels())');
            end
        end
        
        %% Helper functions
        function image = getImage(obj)
            if ~isempty(obj.pixels)
                image = obj.image;
            else
                % Create list of image IDs
                ids = java.util.ArrayList();
                ids.add(java.lang.Long(obj.imageID));
                
                % Create parameters
                param = omero.sys.ParametersI();
                param.acquisitionData;
                
                % Retrieve list of images and select first one
                list = obj.getSession().getContainerService().getImages('omero.model.Image', ids, param);
                image = list.get(0);
            end
        end
        
        function pixels = getPixels(obj)
            if ~isempty(obj.pixels)
                pixels = obj.pixels;
            else
                % Retrieve pixels ID
                pixelsId = obj.getImage().getPixels(0).getId.getValue;
                
                % Get PixelsI object
                pixelsService=obj.getSession().getPixelsService();
                pixels=pixelsService.retrievePixDescription(pixelsId);
            end
        end
        
        function delete(obj)
            if ~isempty(obj.session)
                obj.session.close()
            end
        end
    end
end