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
            
            % Set the colormap
            colormap(hAxes,obj.Colormap);
            
            % Set the colorbar
            hCbar = findobj(get(hAxes,'Parent'),'Tag','Colorbar');
            if strcmp(obj.Colorbar,'on')
                axis image
                if isempty(hCbar)
                    set(hAxes,'Position',[0.05 0.05 .9 .9]);   
                    hCbar = colorbar('peer',hAxes,'FontSize',12);
                    ylabel(hCbar,obj.Units,'FontSize',12);
                else
                    ylabel(hCbar,obj.Units,'FontSize',12);
                end
            else
                if ~isempty(hCbar),colorbar(hCbar,'delete'); end
                set(hAxes,'XLim',[0 size(data,2)],'YLim',[0 size(data,1)],...
                'Position',[0 0 1 1]);
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