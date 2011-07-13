classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        Colormap ='gray';
        Colorbar ='off';
        CLim = [];
        Units='';
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
            % Plot the image and associate the tag
            h=imshow(data,varargin{:});
            set(h,'Tag',tag);
            
            % Clean existing image and set image at the bottom of the stack
            hAxes = get(h,'Parent');
            child=get(hAxes,'Children');
            imChild = child(strcmp(get(child,'Type'),'image'));
            delete(imChild(imChild~=h));
            uistack(h,'bottom');
            
            axis manual
            % Set the colormap
            colormap(hAxes,obj.Colormap);
            
            % Set the colorbar
            hCbar = findobj(get(hAxes,'Parent'),'Tag','Colorbar');
            if strcmp(obj.Colorbar,'on') && isempty(hCbar)
                set(gca,'CLimMode','auto');
                hCBar = colorbar('peer',hAxes);
                ylabel(hCBar,obj.Units);
            elseif strcmp(obj.Colorbar,'off') && ~isempty(hCbar)
                colorbar(hCbar,'delete');
            end
            
            % Set the color limits
            if ~isempty(obj.CLim),set(hAxes,'CLim',obj.CLim); end
        end
        function updateDraw(obj,h,data)
            set(h,'CData',data)
        end
        
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Colormap',obj.Colormap,@ischar);
            ip.addParamValue('Colorbar',obj.Colorbar,@ischar);
            ip.addParamValue('CLim',obj.CLim,@(x) isempty(x) ||isvector(x));
        end
        
       function setProperties(obj,ip)
            obj.Colormap=ip.Results.Colormap;
            obj.Colorbar=ip.Results.Colorbar;
            obj.CLim=ip.Results.CLim;
        end
        
    end 
   
 
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end