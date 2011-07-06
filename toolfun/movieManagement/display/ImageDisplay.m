classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        color = [1 1 1];
        colormap
    end
    methods
        function h=initDraw(obj,data,tag,varargin)
            h=imshow(data,varargin{:});
            set(h,'Tag',tag);
            
            % Clean existing image and set image at the bottom of the stack
            child=get(get(h,'Parent'),'Children');
            imChild = child(strcmp(get(child,'Type'),'image'));
            delete(imChild(imChild~=h));
            uistack(h,'bottom');
        end
        function updateDraw(obj,h,data)
            set(h,'CData',data)
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.color,@isvector);
        end
        
       function setProperties(obj,ip)
            obj.color=ip.Results.Color;
        end
        
    end 
   
 
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end