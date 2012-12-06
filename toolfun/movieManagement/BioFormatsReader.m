classdef  BioFormatsReader < Reader
    % Concrete implementation of MovieObject for a single movie
    
    properties (Transient =true)
        formatReader
    end
    
    methods
        %% Constructor
        function obj = BioFormatsReader(path, iSeries)
            bfCheckJavaPath(); % Check loci-tools.jar is in the Java path
            loci.common.DebugTools.enableLogging('OFF');
            obj.formatReader = bfGetReader(path, false);
            if nargin>1,
                obj.formatReader.setSeries(iSeries);
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
        
        function fileNames = getImageFileNames(obj, iChan, varargin)
            % Generate image file names
            [~, fileName] = fileparts(char(obj.formatReader.getCurrentFile));
            basename = sprintf('%s_s%g_c%d_t',fileName, obj.getSeries()+1, iChan);
            fileNames = arrayfun(@(t) [basename num2str(t, ['%0' num2str(floor(log10(obj.getSizeT))+1) '.f']) '.tif'],...
                1:obj.getSizeT,'Unif',false);
        end
        
        function chanNames = getChannelNames(obj, iChan)
            [~, fileName, fileExt] = fileparts(char(obj.formatReader.getCurrentFile));
            channelID = @(x) char(obj.getMetadataStore().getChannelID(obj.getSeries(), x-1));
            chanNames = arrayfun(@(x) [fileName fileExt ':'  channelID(x)], iChan, ...
                'Unif',false);
        end
        
        function I = loadImage(obj, c, t)
            % Using bioformat tools, get the reader and retrieve dimension order
            r = obj.formatReader;
            I = zeros([obj.getSizeY(), obj.getSizeX(), numel(t)]);
            z = 1;
            for i = 1 : numel(t),
                iPlane = loci.formats.FormatTools.getIndex(r, z-1, c-1, t(i)-1);
                I(:,:,i) = bfGetPlane(r, iPlane + 1);
            end
        end
        
        function delete(obj)
            obj.formatReader.close()
        end
    end
end