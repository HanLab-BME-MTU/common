classdef  BioFormatsReader < Reader
    % Concrete implementation of MovieObject for a single movie
    
    properties (Transient =true)
        formatReader
    end
    
    methods
        %% Constructor
        function obj = BioFormatsReader(varargin)
            % Check loci-tools.jar is in the Java path
            if isa(varargin{1}, 'loci.formats.IFormatReader'),
               obj.formatReader = varargin{1};
            else
                loci.common.DebugTools.enableLogging('OFF');
                obj.formatReader = bfGetReader(varargin{1}, false);
            end
            if nargin>1
                obj.formatReader.setSeries(varargin{2});
            end
        end
        
        function metadataStore = getMetadataStore(obj)
            r = obj.formatReader;
            metadataStore = r.getMetadataStore();
        end
        
        function series = getSeries(obj)
            series = obj.formatReader.getSeries();
        end
        
        function sizeX = getSizeX(obj, varargin)
            sizeX = obj.getMetadataStore().getPixelsSizeX(obj.getSeries()).getValue();
        end
        
        function sizeY = getSizeY(obj, varargin)
            sizeY = obj.getMetadataStore().getPixelsSizeY(obj.getSeries()).getValue();
        end
        
        function sizeZ = getSizeZ(obj, varargin)
            sizeZ = obj.getMetadataStore().getPixelsSizeZ(obj.getSeries()).getValue();
        end
        
        function sizeT = getSizeT(obj, varargin)
            sizeT = obj.getMetadataStore().getPixelsSizeT(obj.getSeries()).getValue();
        end
        
        function sizeC = getSizeC(obj, varargin)
            sizeC = obj.getMetadataStore().getPixelsSizeC(obj.getSeries()).getValue();
        end
        
        function bitDepth = getBitDepth(obj, varargin)
            pixelType = obj.formatReader.getPixelType();
            bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
            bitDepth = 8 * bpp;
        end
        
        function fileNames = getImageFileNames(obj, iChan, varargin)
            % Generate image file names
            [~, fileName] = fileparts(char(obj.formatReader.getCurrentFile));
            basename = sprintf('%s_s%g_c%d_t',fileName, obj.getSeries()+1, iChan);
            fileNames = arrayfun(@(t) [basename num2str(t, ['%0' num2str(floor(log10(obj.getSizeT))+1) '.f']) '.tif'],...
                1:obj.getSizeT,'Unif',false);
        end
        
        function channelNames = getChannelNames(obj, iChan)
            [~, fileName, fileExt] = fileparts(char(obj.formatReader.getCurrentFile));
            
            if obj.formatReader.getSeriesCount() > 1
                base = [fileName fileExt ' Series ' num2str(obj.getSeries()+1) ' Channel '];
            else
                base = [fileName fileExt ' Channel '];
            end
            
            channelNames = arrayfun(@(x) [base num2str(x)], iChan, 'Unif',false);
        end
        
        function index = getIndex(obj, z, c, t)
            index = loci.formats.FormatTools.getIndex(obj.formatReader, z, c, t);
        end
        
        function I = loadImage(obj, c, t, varargin)
            
            ip = inputParser;
            ip.addRequired('c', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC()));
            ip.addRequired('t', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT()));
            ip.addOptional('z', 1, @(x) isscalar(x) && ismember(x, 1 : obj.getSizeZ()));
            ip.parse(c, t, varargin{:});
            
            % Using bioformat tools, get the reader and retrieve dimension order
            javaIndex =  obj.getIndex(ip.Results.z - 1, c - 1, t - 1);
            I = bfGetPlane(obj.formatReader, javaIndex + 1);
        end
        
        function delete(obj)
            obj.formatReader.close()
        end
    end
end