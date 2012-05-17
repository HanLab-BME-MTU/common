classdef ImageDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        Colormap ='gray';
        Colorbar ='off';
        ColorbarLocation ='EastOutside';
        CLim = [];
        Units='';
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
        ScaleFactor = 1;
    end
    methods
        function obj=ImageDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
            
        function h=initDraw(obj,data,tag,varargin)
            % Plot the image and associate the tag
            h=imshow(data/obj.ScaleFactor,varargin{:});
            set(h,'Tag',tag,'CDataMapping','scaled');
            hAxes = get(h,'Parent');
            set(hAxes,'XLim',[0 size(data,2)],'YLim',[0 size(data,1)]);
            obj.applyImageOptions(h)
        end
        function updateDraw(obj,h,data)
            set(h,'CData',data/obj.ScaleFactor)
            obj.applyImageOptions(h)
        end
        
        function applyImageOptions(obj,h)
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
            axesPosition = [0 0 1 1];
            if strcmp(obj.Colorbar,'on')
                if length(obj.ColorbarLocation) >6 && strcmp(obj.ColorbarLocation(end-6:end),'Outside'),
                    axesPosition = [0.05 0.05 .9 .9];
                end
                if isempty(hCbar)
                    set(hAxes,'Position',axesPosition);   
                    hCbar = colorbar('peer',hAxes,obj.sfont{:});
                end
                set(hCbar,'Location',obj.ColorbarLocation);
                ylabel(hCbar,obj.Units,obj.lfont{:});
            else
                if ~isempty(hCbar),colorbar(hCbar,'delete'); end
                set(hAxes,'Position',axesPosition);
            end
            
            % Set the color limits
            if ~isempty(obj.CLim),set(hAxes,'CLim',obj.CLim/obj.ScaleFactor); end
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
            params(7).name='ColorbarLocation';
            findclass(findpackage('scribe'),'colorbar');
            locations = findtype('ColorbarLocationPreset');
            locations = locations.Strings;
            params(7).validator=@(x) any(strcmp(x,locations));
            params(8).name='ScaleFactor';
            params(8).validator=@isscalar;
        end
        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end