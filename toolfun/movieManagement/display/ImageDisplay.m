classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        color = '';
        colormap
    end
    methods
        function h=initDraw(obj,data,tag,varargin)
            rgbData=repmat(double(data)/max(data(:)),[1 1 3]);
            switch obj.color
                case 'r'
                    rgbData(:,:,2:3)=0;
                case 'b'
                    rgbData(:,:,1:2)=0;
                case 'g'
                    rgbData(:,:,[1 3])=0;
            end
            h=imshow(rgbData,varargin{:});
            set(h,'Tag',tag);
        end
        function setProperties(obj,ip)
            obj.color=ip.Results.Color;
        end
        function updateDraw(obj,h,data)
            rgbData=get(h,'CData');
            switch obj.color
                case 'r'
                    rgbData(:,:,1)=double(data)/max(data(:));
                case 'g'
                    rgbData(:,:,2)=double(data)/max(data(:));
                case 'b'
                    rgbData(:,:,3)=double(data)/max(data(:));
                otherwise
                    rgbData=repmat(double(data)/max(data(:)),[1 1 3]);
            end
            set(h,'CData',rgbData)
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.color,@ischar);
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end