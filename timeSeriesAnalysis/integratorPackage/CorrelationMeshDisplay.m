classdef CorrelationMeshDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        LineStyle = '-'
        Color='r';        
    end
    properties (SetAccess = protected)
        slider;
    end
    methods
        function obj=MeshDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
  
            delete(gca);
            [nx,ny,nz,nz2] = size(data.Z);
            if ~isfield(data,'X'), data.X=repmat(1:nx,1,ny); end
            if ~isfield(data,'Y'), data.Y= repmat(1:ny,nx,1); end
            h(1)=mesh(data.X(:,:,1),data.Y,data.Z(:,:,1),'FaceColor','interp');
           
            xlabel('Lag (s)')
            ylabel('Window number')
            zlabel('Correlation')
            
            if ~isempty(data.bounds)
                hold on
%                 [cL,~] = size(data.Z);
                upline  = repmat(data.bounds(1,:,1),nx,1);
                h(2)=mesh(data.X(:,:,1),data.Y,upline,'FaceColor',[63/255 63/255 63/255]);
                
                dline  = repmat(data.bounds(2,:,1),nx,1);
                h(3)=mesh(data.X(:,:,1),data.Y,dline,'FaceColor',[63/255 63/255 63/255]);
            end
            
            if nz>1
                axesPos = get(get(h(1),'Parent'),'Position');
                mainFig = get(get(h(1),'Parent'),'Parent');
                set(mainFig,'Toolbar','figure');
                obj.slider = uicontrol(mainFig,'Style','slider',...
                    'Units','normalized',...
                    'Position',[axesPos(1) 1-.03 axesPos(3) .03],...
                    'Value',1,'Min',1,'Max',nz,...
                    'SliderStep',[1/(nz-1)  5/(nz-1)],...
                    'Tag','slider_depth','BackgroundColor','white',...
                    'Callback',@(hObject,event) updateDraw(obj,h,data));
            end
            
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            [nx,ny,nz] = size(data.Z);
            if nz>1,
                depth = round(get(obj.slider,'Value'));
            else
                depth=1;
            end
            set(h(1),'XData',data.X(:,:,depth),'ZData',data.Z(:,:,depth));
            if ~isempty(data.bounds)
                upline  = repmat(data.bounds(1,:,depth),nx,1);
                set(h(2),'XData',data.X(:,:,depth),'ZData',upline);                
                dline  = repmat(data.bounds(2,:,depth),nx,1);
                set(h(3),'XData',data.X(:,:,depth),'ZData',dline);
            end

        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@ischar);
            ip.addParamValue('Marker',obj.Marker,@ischar);
            ip.addParamValue('LineStyle',obj.LineStyle,@ischar);  
        end 
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.Marker=ip.Results.Marker;
            obj.LineStyle=ip.Results.LineStyle;
        end
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isstruct;
        end
    end    
end