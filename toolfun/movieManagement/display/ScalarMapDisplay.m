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
            
            if size(data,3)>1    
                dataSlice = squeeze(data(:,1,:));
                h=imagesc(dataSlice,varargin{:});
                mainFig = get(get(h,'Parent'),'Parent');
                obj.slider = uicontrol(mainFig,'Style','slider',...
                    'Position',[20 20 30 250],...
                    'Value',1,'Min',1,'Max',size(data,2),...
                    'SliderStep',[1/(size(data,2)-1)  5/(size(data,2)-1)],...
                    'Tag','slider_depth','BackgroundColor','white',...
                    'Callback',@(hObject,event) updateDraw(obj,h,data));
                set(h,'AlphaData',~isnan(dataSlice))
            else
                h=imagesc(data,varargin{:});
                set(h,'alphadata',~isnan(data))
            end

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
            if ~isempty(obj.Labels{1}),xlabel(obj.Labels{1}); end
            if ~isempty(obj.Labels{2}),ylabel(obj.Labels{2}); end
        end

        function updateDraw(obj,h,data)
            if size(data,3)>1
                depth = round(get(obj.slider,'Value'));
                dataSlice=squeeze(data(:,depth,:));
                set(h,'CData',dataSlice);
                set(h,'AlphaData',~isnan(dataSlice));
            else
                set(h,'CData',data);
            end
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