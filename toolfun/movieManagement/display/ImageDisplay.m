classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        Colormap ='gray';
        Colorbar ='off';
        CLim = [];
        Units='';
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
    end
    methods
        function obj=ImageDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
            
        function h=initDraw(obj,data,tag,varargin)
            % Plot the image and associate the tag
            h=imshow(data,varargin{:});
            set(h,'Tag',tag,'CDataMapping','scaled');
            hAxes = get(h,'Parent');
            set(hAxes,'XLim',[0 size(data,2)],'YLim',[0 size(data,1)]);
            obj.applyImageOptions(h,data)
        end
        function updateDraw(obj,h,data)
            set(h,'CData',data)
            obj.applyImageOptions(h,data)
        end
        
        function applyImageOptions(obj,h,data)
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
                if isempty(hCbar)
                    set(hAxes,'Position',[0.05 0.05 .9 .9]);   
                    hCbar = colorbar('peer',hAxes,obj.sfont{:});
                end
                ylabel(hCbar,obj.Units,obj.lfont{:});
            else
                if ~isempty(hCbar),colorbar(hCbar,'delete'); end
                set(hAxes,'Position',[0 0 1 1]);
            end

            
            % Set the color limits
            if ~isempty(obj.CLim),set(hAxes,'CLim',obj.CLim); end
        end
    end 
 
    methods (Static)
         function params=getParamValidators()
            params(1).name='Colormap';
            params(1).validator=@ischar;
            params(2).name='Colorbar';
            params(2).validator=@(x) any(strcmp(x,{'on','off'}));
            params(3).name='CLim';
            params(3).validator=@isvector;
            params(4).name='Units';
            params(4).validator=@ischar;
            params(5).name='sfont';
            params(5).validator=@iscell;
            params(6).name='lfont';
            params(6).validator=@iscell;
        end
        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end