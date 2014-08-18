classdef  BioFormatsReader < Reader
    % Concrete implementation of MovieObject for a single movie
    
    properties (Transient =true)
        formatReader
        series = 0;
    end
    
    methods
        %% Constructor
        function obj = BioFormatsReader(varargin)
            % Check loci-tools.jar is in the Java path
            bfCheckJavaPath();
            if isa(varargin{1}, 'loci.formats.IFormatReader'),
                obj.formatReader = varargin{1};
                obj.series = obj.formatReader.getSeries();
            else
                loci.common.DebugTools.enableLogging('OFF');
                obj.formatReader = bfGetReader(varargin{1}, false);
            end
            if nargin>1,
                obj.series = varargin{2};
                obj.formatReader.setSeries(obj.series);
            end
        end
        
        function metadataStore = getMetadataStore(obj)
            metadataStore = obj.formatReader.getMetadataStore();
        end
        
        function r = getReader(obj)
            r = obj.formatReader;
            r.setSeries(obj.getSeries());
        end
        
        function series = getSeries(obj)
            series = obj.series;
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
            pixelType = obj.getReader().getPixelType();
            bpp = loci.formats.FormatTools.getBytesPerPixel(pixelType);
            bitDepth = 8 * bpp;
        end
        
        function fileNames = getImageFileNames(obj, iChan, varargin)
            % Generate image file names
            usedFiles = obj.getReader().getUsedFiles();
            [~, fileName] = fileparts(char(usedFiles(1)));
            basename = sprintf('%s_s%g_c%d_t',fileName, obj.getSeries()+1, iChan);
            fileNames = arrayfun(@(t) [basename num2str(t, ['%0' num2str(floor(log10(obj.getSizeT))+1) '.f']) '.tif'],...
                1:obj.getSizeT,'Unif',false);
        end
        
        function channelNames = getChannelNames(obj, iChan)
            usedFiles = obj.getReader().getUsedFiles();
            [~, fileName, fileExt] = fileparts(char(usedFiles(1)));
            base = [fileName fileExt];
            if obj.getReader().getSeriesCount() > 1
                base = [base ' Series ' num2str(obj.getSeries()+1)];
            end
            base = [base ' Channel '];
            
            channelNames = arrayfun(@(x) [base num2str(x)], iChan, 'Unif',false);
        end
        
        function index = getIndex(obj, z, c, t)
            index = loci.formats.FormatTools.getIndex(obj.getReader(), z, c, t);
        end
        
        function I = loadImage(obj, c, t, varargin)
            % Retrieve single plane specified by its (c, t, z) coordinates
            
            ip = inputParser;
            ip.addRequired('c', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC()));
            ip.addRequired('t', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT()));
            ip.addOptional('z', 1, @(x) isscalar(x) && ismember(x, 1 : obj.getSizeZ()));
            ip.parse(c, t, varargin{:});
            
            % Using bioformat tools, get the reader and retrieve dimension order
            javaIndex =  obj.getIndex(ip.Results.z - 1, c - 1, t - 1);
            I = bfGetPlane(obj.getReader(), javaIndex + 1);
        end
        
        function I = loadStack(obj, c, t, varargin)
            % Retrieve entire z-stack or sub-stack
            
            % Input check
            ip = inputParser;
            ip.addRequired('c', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeC()));
            ip.addRequired('t', @(x) isscalar(x) && ismember(x, 1 : obj.getSizeT()));
            ip.addOptional('z', 1 : obj.getSizeZ(), @(x) all(ismember(x, 1 : obj.getSizeZ())));
            ip.parse(c, t, varargin{:});
            
            % Determine image class from pixel type
            r = obj.getReader();
            pixelType = r.getPixelType();
            pixelTypeString = loci.formats.FormatTools.getPixelTypeString(pixelType);
            if strcmp(char(pixelTypeString), 'float'),
                % Handle float/single conversion
                imClass = 'single';
            else
                imClass = char(pixelTypeString);
            end
            % Load image stack by looping over bfGetPlane
            I = zeros(obj.getSizeY(), obj.getSizeX(),numel(ip.Results.z), imClass);
            for iz = 1 : numel(ip.Results.z)
                javaIndex =  obj.getIndex(ip.Results.z(iz) - 1, c - 1, t - 1);
                I(:, :, iz) = bfGetPlane(obj.getReader(), javaIndex + 1);
            end
        end
        
        function delete(obj)
            obj.formatReader.close()
        end
    end
end