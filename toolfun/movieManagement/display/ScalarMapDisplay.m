classdef ScalarMapDisplay < MovieDataDisplay
    %Abstract class for displaying image processing output
    properties
        Colormap='jet';
        Colorbar ='on';
        CLim = [];
        ScaleLabel='';
        Labels={'',''};
        depthDim=3;
    end
    properties (SetAccess = protected)
        slider;
    end
    
    
    methods
        function obj=ScalarMapDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        
        function h=initDraw(obj,data,tag,varargin)
            
            h=imagesc(data(:,:,1),varargin{:});
            set(h,'AlphaData',~isnan(data(:,:,1)))
            % Plot the image and associate the tag

            set(h,'Tag',tag,'UserData',data);
            
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
%                     set(hAxes,'Position',[0.05 0.05 .9 .9]);   
                    hCBar = colorbar('peer',hAxes,'FontSize',12);
                    ylabel(hCBar,obj.ScaleLabel,'FontSize',12);
                end
            else
                if ~isempty(hCbar),colorbar(hCbar,'delete'); end
                set(hAxes,'XLim',[0 size(data,2)],'YLim',[0 size(data,1)],...
                'Position',[0 0 1 1]);
            end
            
            % Set the color limits
            if ~isempty(obj.CLim),set(hAxes,'CLim',obj.CLim); end
                        
            % Set the color limits
            if ~isempty(obj.Labels{1}),xlabel(obj.Labels{1},'Parent',hAxes); end
            if size(data,3)>1   
                if ~isempty(obj.Labels{3}),ylabel(obj.Labels{3},'Parent',hAxes); end
            else
                if ~isempty(obj.Labels{2}),ylabel(obj.Labels{2},'Parent',hAxes); end
            end
            
            nz=size(data,3);
            if nz>1
                axesPos = get(hAxes,'Position');
                mainFig = get(get(h,'Parent'),'Parent');
                set(mainFig,'Toolbar','figure');
                obj.slider = uicontrol(mainFig,'Style','slider',...
                    'Units','normalized',...
                    'Position',[axesPos(1)-.05 axesPos(2) axesPos(3)/20 axesPos(4)],...
                    'Value',1,'Min',1,'Max',nz,'SliderStep',[1/(nz-1)  5/(nz-1)],...
                    'Tag','slider_depth','BackgroundColor','white',...
                    'Callback',@(hObject,event) updateDraw(obj,h,data));
                if ~isempty(obj.Labels{2}),
                    hp = uipanel(mainFig, 'Units','normalized',...
                    'Position',[axesPos(1)/4 axesPos(2) axesPos(3)/20 axesPos(4)],...
                    'BorderType','none',...
                    'BackgroundColor',get(gcf,'Color')); 
                    ha = axes('Parent',hp,'Visible','off');
                    text(0,1/3,obj.Labels{2},'Parent',ha,'rotation',90)
                end
            end
        end

        function updateDraw(obj,h,data)
            if size(data,3)>1
                depth = round(get(obj.slider,'Value'));
            else
                depth=1;
            end
            set(h,'CData',data(:,:,depth));
            set(h,'AlphaData',~isnan(data(:,:,depth)));
            
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