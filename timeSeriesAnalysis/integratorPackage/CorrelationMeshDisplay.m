classdef CorrelationMeshDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        LineStyle = '-'
        Color='r';        
    end
    properties (SetAccess = protected)
        slider;
        slider2;
    end
    methods
        function obj=MeshDisplay(varargin)
            obj@MovieDataDisplay(varargin{:})
        end
        function h=initDraw(obj,data,tag,varargin)
  
            delete(gca);
            [nx,ny,nz,nz2] = size(data.Z);
            if ~isfield(data,'X'), data.X=repmat(1:nx,1,ny); end
            if ~isfield(data,'Y'), data.Y= repmat(1:ny,nx,1); end
            h(1)=mesh(data.X(:,:,1,1),data.Y,data.Z(:,:,1,1),'FaceColor','interp');
           
            xlabel('Lag (s)')
            ylabel('Window number')
            zlabel('Correlation')
            
            if ~isempty(data.bounds)
                hold on
%                 [cL,~] = size(data.Z);
                upline  = repmat(data.bounds(1,:,1,1),nx,1);
                h(2)=mesh(data.X(:,:,1,1),data.Y,upline,'FaceColor',[63/255 63/255 63/255]);
                
                dline  = repmat(data.bounds(2,:,1,1),nx,1);
                h(3)=mesh(data.X(:,:,1,1),data.Y,dline,'FaceColor',[63/255 63/255 63/255]);
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
            if nz2>1
                axesPos = get(get(h(1),'Parent'),'Position');
                mainFig = get(get(h(1),'Parent'),'Parent');
                set(mainFig,'Toolbar','figure');
                obj.slider2 = uicontrol(mainFig,'Style','slider',...
                    'Units','normalized',...
                    'Position',[axesPos(1) 1-.06 axesPos(3) .03],...
                    'Value',1,'Min',1,'Max',nz,...
                    'SliderStep',[1/(nz2-1)  5/(nz2-1)],...
                    'Tag','slider_depth','BackgroundColor','white',...
                    'Callback',@(hObject,event) updateDraw(obj,h,data));
            end
            
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            [nx,ny,nz,nz2] = size(data.Z);
            if nz>1,
                zdepth = round(get(obj.slider,'Value'));
            else
                zdepth=1;
            end
            if nz2>1,
                zdepth2 = round(get(obj.slider2,'Value'));
            else
                zdepth2=1;
            end
            set(h(1),'XData',data.X(:,:,zdepth,zdepth2),'ZData',data.Z(:,:,zdepth,zdepth2));
            if ~isempty(data.bounds)
                upline  = repmat(data.bounds(1,:,zdepth,zdepth2),nx,1);
                set(h(2),'XData',data.X(:,:,zdepth,zdepth2),'ZData',upline);                
                dline  = repmat(data.bounds(2,:,zdepth,zdepth2),nx,1);
                set(h(3),'XData',data.X(:,:,zdepth,zdepth2),'ZData',dline);
            end
        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@ischar;
            params(2).name='Marker';
            params(2).validator=@ischar;
            params(3).name='Marker';
            params(3).validator=@ischar;
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end    
end