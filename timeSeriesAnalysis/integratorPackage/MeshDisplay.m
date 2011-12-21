classdef MeshDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        LineStyle = '-'
        Color='r';
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
        XLabel='';
        YLabel='';
        ZLabel='';
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
            delete(gca)
            % Read sizes and replicate data if applicable
            [nx,ny,nBands,nBands2] = size(data.Z);
            if isvector(data.X), data.X=repmat(data.X,1,ny); end
            if ~isfield(data,'Y'), data.Y= repmat(1:ny,nx,1); end
            
            % Create mesh display
            h(1)=mesh(data.X,data.Y,data.Z(:,:,1,1),'FaceColor','interp');
            if isfield(data,'bounds') && ~isempty(data.bounds)
                hold on
                upline  = repmat(data.bounds(1,:,1,1),nx,1);
                h(2)=mesh(data.X,data.Y,upline,'FaceColor',[63/255 63/255 63/255]);
                
                dline  = repmat(data.bounds(2,:,1,1),nx,1);
                h(3)=mesh(data.X,data.Y,dline,'FaceColor',[63/255 63/255 63/255]);
                zLim =[min(vertcat(data.Z(:),data.bounds(:))) max(vertcat(data.Z(:),data.bounds(:)))];
            else
                zLim =[min(data.Z(:)) max(data.Z(:))];
            end
            set(h,'Tag',tag);
            
            % Set axis options
            xLim=[min(data.X(:)) max(data.X(:))];
            yLim=[min(data.Y(:)) max(data.Y(:))];
            xlabel(obj.XLabel,obj.lfont{:});
            ylabel(obj.YLabel,obj.lfont{:});
            zlabel(obj.ZLabel,obj.lfont{:})
            set(gca,'Linewidth',1.5,obj.sfont{:},'XLim',xLim,'YLim',yLim,'ZLim',zLim);
            
            % Create sliders if multiple bands
            axesPos = get(get(h(1),'Parent'),'Position');
            mainFig = get(get(h(1),'Parent'),'Parent');

            if nBands>1
                set(mainFig,'Toolbar','figure');
                obj.slider = uicontrol(mainFig,'Style','slider',...
                    'Units','normalized',...
                    'Position',[axesPos(1) 1-.03 axesPos(3) .03],...
                    'Value',1,'Min',1,'Max',nBands,...
                    'SliderStep',[1/(nBands-1)  5/(nBands-1)],...
                    'BackgroundColor','white',...
                    'Callback',@(hObject,event) updateDraw(obj,h,data));
            end
            if nBands2>1
                set(mainFig,'Toolbar','figure');
                obj.slider2 = uicontrol(mainFig,'Style','slider',...
                    'Units','normalized',...
                    'Position',[axesPos(1) 1-.06 axesPos(3) .03],...
                    'Value',1,'Min',1,'Max',nBands2,...
                    'SliderStep',[1/(nBands2-1)  5/(nBands2-1)],...
                    'BackgroundColor','white',...
                    'Callback',@(hObject,event) updateDraw(obj,h,data));
            end

            
        end
        function updateDraw(obj,h,data)
            [nx,ny,nBands,nBands2] = size(data.Z);
            if nBands>1,
                zdepth = round(get(obj.slider,'Value'));
            else
                zdepth=1;
            end
            if nBands2>1,
                zdepth2 = round(get(obj.slider2,'Value'));
            else
                zdepth2=1;
            end
            set(h(1),'XData',data.X,'ZData',data.Z(:,:,zdepth,zdepth2));
            if isfield(data,'bounds') && ~isempty(data.bounds)
                upline  = repmat(data.bounds(1,:,zdepth,zdepth2),nx,1);
                set(h(2),'XData',data.X,'ZData',upline);                
                dline  = repmat(data.bounds(2,:,zdepth,zdepth2),nx,1);
                set(h(3),'XData',data.X,'ZData',dline);
            end
        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@ischar;
            params(2).name='Marker';
            params(2).validator=@ischar;
            params(3).name='LineStyle';
            params(3).validator=@ischar;
            params(4).name='sfont';
            params(4).validator=@iscell;
            params(5).name='lfont';
            params(5).validator=@iscell;
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end    
end