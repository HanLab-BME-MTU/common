classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        color = [1 1 1];
        colormap
    end
    methods
        function h=initDraw(obj,data,tag,varargin)
            h=imshow(obj.formatData(data),varargin{:});
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            set(h,'CData',obj.formatData(data))
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.color,@isvector);
        end
        
       function setProperties(obj,ip)
            obj.color=ip.Results.Color;
        end
        
    end 
    
    methods(Access=protected)
        function rgbData=formatData(obj,data)
            data =double(data);
            rgbData=repmat(data/max(data(:)),[1 1 3]);
            for  i=1:3
                rgbData(:,:,i)=rgbData(:,:,i)*obj.color(i);
            end
        end
    end
 
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end