classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        color = [1 1 1];
        Colormap ='gray';
    end
    methods
        function obj=ImageDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            h=imshow(data,varargin{:});
            set(h,'Tag',tag);
            
            % Clean existing image and set image at the bottom of the stack
            child=get(get(h,'Parent'),'Children');
            imChild = child(strcmp(get(child,'Type'),'image'));
            delete(imChild(imChild~=h));
            uistack(h,'bottom');
            
            % Set the colormap
            colormap(get(h,'Parent'),obj.Colormap);
%             set(get(get(h,'Parent'),'Parent'),'Colormap',obj.colormap);
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