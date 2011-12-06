classdef ScalarMapDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        Colormap='jet';
        NaNColor = [1 1 1];
        Colorbar ='on';
        CLim = [];
        Units='';
        Labels={'',''};
        depthDim=3;
    end
    properties (SetAccess = protected)
        slider;
    end

    methods
        function obj=ScalarMapDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
        
        function h=initDraw(obj,data,tag,varargin)
            % Create extended cmap (for NaNs)
            cmap = [obj.NaNColor;colormap(obj.Colormap)];
            imData= data(:,:,1);
            imData(isnan(imData)) = min(data(:))-(max(data(:)-min(data(:))))*1e-10;
            h=imagesc(imData,varargin{:});
            
            % Plot the image and associate the tag
            set(h,'Tag',tag,'UserData',data);
            
            % Clean existing image and set image at the bottom of the stack
            hAxes = get(h,'Parent');
            child=get(hAxes,'Children');
            imChild = child(strcmp(get(child,'Type'),'image'));
            delete(imChild(imChild~=h));
            uistack(h,'bottom');
            
            % Set the colormap
            colormap(hAxes,cmap);
            
            % Set the fonts
            sfont = {'FontName', 'Helvetica', 'FontSize', 18};
            lfont = {'FontName', 'Helvetica', 'FontSize', 22};
            
            % Set the colorbar
            hCbar = findobj(get(hAxes,'Parent'),'Tag','Colorbar');
            if strcmp(obj.Colorbar,'on')
                axis image
                if isempty(hCbar)
                    %  set(hAxes,'Position',[0.05 0.05 .9 .9]);
                    hCBar = colorbar('peer',hAxes,'FontSize',12);
                    ylabel(hCBar,obj.Units,lfont{:});
                end
            else
                if ~isempty(hCbar),colorbar(hCbar,'delete'); end
                set(hAxes,'XLim',[0 size(data,2)],'YLim',[0 size(data,1)],...
                'Position',[0 0 1 1]);
            end
            
            % Set the color limits
            if ~isempty(obj.CLim),set(hAxes,'CLim',obj.CLim); end
                        
            % Set the color limits
            if ~isempty(obj.Labels{1}),xlabel(obj.Labels{1},'Parent',hAxes,lfont{:}); end
            if size(data,3)>1   
                if ~isempty(obj.Labels{3}),ylabel(obj.Labels{3},'Parent',hAxes,lfont{:}); end
            else
                if ~isempty(obj.Labels{2}),ylabel(obj.Labels{2},'Parent',hAxes,lfont{:}); end
            end
            set(hAxes,'LineWidth', 1.5, sfont{:})
            nz=size(data,3);
            if nz>1

                axesPos = get(hAxes,'Position');
                mainFig = get(get(h,'Parent'),'Parent');
                set(mainFig,'Toolbar','figure');
                obj.slider = uicontrol(mainFig,'Style','slider',...
                    'Units','normalized',...
                    'Position',[axesPos(1) 1-.06 axesPos(3) .03],...
                    'Value',1,'Min',1,'Max',nz,'SliderStep',[1/(nz-1)  5/(nz-1)],...
                    'Tag','slider_depth','BackgroundColor','white',...
                    'Callback',@(hObject,event) updateDraw(obj,h,data));
%                     'Position',[axesPos(1)-.05 axesPos(2) axesPos(3)/20
%                     axesPos(4)],...
                
                if ~isempty(obj.Labels{2}),
                    uicontrol(mainFig, 'Style','text', 'Units','normalized',...
                    'Position',[axesPos(1) 1-.03 axesPos(3) .03],...
                    'String',obj.Labels{2},...
                    'BackgroundColor',get(gcf,'Color')); 
                end
            end
        end

        function updateDraw(obj,h,data)
            if size(data,3)>1
                depth = round(get(obj.slider,'Value'));
            else
                depth=1;
            end
            
            imData= data(:,:,depth);
            imData(isnan(imData)) = 0;

            set(h,'CData',imData);      
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
            params(5).name='NaNColor';
            params(5).validator=@(x) isequal(size(x),[1 3]);
        end

        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end