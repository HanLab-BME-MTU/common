classdef MockReader < Reader
    properties
        images = {}
    end
    methods
        function obj = MockReader(varargin)
            if(nargin < 1)
                obj.images = obj.createDefaultData();
            else
                obj.images = varargin{1};
            end
        end
        function I = createDefaultData(obj)
            sizeY = 10;
            sizeX = 20;
            sizeC = 3;
            sizeT = 5;
            sizeZ = 7;
            % create channels
            I = cell(1,sizeC);
            for c = 1 : sizeC
                I{c} = cell(1,sizeT);
                for t = 1 : sizeT
                    I{c}{t} = cell(1,sizeZ);
                    for z = 1 : sizeZ
                        I{c}{t}{z} = ones(sizeY,sizeX,'uint16');
                        I{c}{t}{z}(:) = [ c t z]*[100 10 1]';
                    end
                end
            end
        end
        function o = getSizeX(obj,c,t,z,varargin)
            if(nargin < 2); c = 1; end;
            if(nargin < 3); t = 1; end;
            if(nargin < 4); z = 1; end;
            o = size(obj.images{c}{t}{z},2);
        end
        function o = getSizeY(obj,c,t,z,varargin)
            if(nargin < 2); c = 1; end;
            if(nargin < 3); t = 1; end;
            if(nargin < 4); z = 1; end;
            o = size(obj.images{c}{t}{z},1);
        end
        function o = getSizeZ(obj,c,t,varargin)
            if(nargin < 2); c = 1; end;
            if(nargin < 3); t = 1; end;
            o = size(obj.images{c}{t},2);
        end
        function o = getSizeC(obj,varargin)
            o = length(obj.images);
        end
        function o = getSizeT(obj,c,varargin)
            if(nargin < 2); c = 1; end;
            o = size(obj.images{c},2);
        end
        function o = getBitDepth(obj,c,t,z,varargin)
            if(nargin < 2); c = 1; end;
            if(nargin < 3); t = 1; end;
            if(nargin < 4); z = 1; end;
            o = str2num( strrep(class(obj.images{c}{t}{z}),'uint','') );
        end
        function o = getImageFileNames(obj,c,t,z,varargin)
            if(nargin < 2); c = 1; end;
            if(nargin < 3); t = 1; end;
            if(nargin < 4); z = 1; end;
            o = ['mockReader' num2str(obj.images{c}{t}{z}(1,1)) '.test'];
        end
        function o = getChannelNames(obj,c)
            if(nargin < 2); c = 1; end;
            o = ['ch' num2str(c) ];
        end
        function I = loadImage(obj,c,t,z)
            if(nargin < 2); c = 1; end;
            if(nargin < 3); t = 1; end;
            if(nargin < 4); z = 1; end;
            I = obj.images{c}{t}{z};
        end
    end
end
