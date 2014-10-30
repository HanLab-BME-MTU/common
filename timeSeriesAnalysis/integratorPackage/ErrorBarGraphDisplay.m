classdef ErrorBarGraphDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        LineStyle = '-';
        LineWidth = 2;
        Color='r';
        XLabel ='';
        YLabel ='';
        Input1='';
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
        nBandMax=5;
    end
    properties (SetAccess = protected)
        bands;
    end
    methods
        function obj=ErrorBarGraphDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
        function h=initDraw(obj,data,tag,varargin)
            
            nBands = min(size(data.Y,2),obj.nBandMax);
            nBands2 = min(size(data.Y,3),obj.nBandMax);
            
            % Plot data as color-coded lines with error bars
            nBandsTot = nBands*nBands2;
            colors = hsv(nBandsTot);
            h=-ones(nBandsTot,3);
            hold on
            for j=1:nBands2
                for i=1:nBands                
                    ind=sub2ind([nBands nBands2],i,j);
                    h(ind,1)=plot(data.X,data.Y(:,i,j),'Linestyle',obj.LineStyle,...
                        'LineWidth',obj.LineWidth,'Color',colors(ind,:));
%                     h(ind,2)=errorbar(data.X,data.Y(:,i,j),...
%                         data.bounds(2,:,i,j)'-data.Y(:,i,j),...
%                         data.bounds(1,:,i,j)'-data.Y(:,i,j),...
%                         'LineWidth', 2,'Color',colors(ind,:));
                    h(ind,2)=plot(data.X,data.bounds(1,:,i,j)',...
                        'LineWidth', 1,'Linestyle','--','Color',colors(ind,:));
                    h(ind,3)=plot(data.X,data.bounds(2,:,i,j)',...
                        'LineWidth', 1,'Linestyle','--','Color',colors(ind,:));
                    
                end
            end
            set(h,'Tag',tag);
            
            % Set axis options
            xlabel(obj.XLabel,obj.lfont{:});
            ylabel(obj.YLabel,obj.lfont{:});
            xLim=[min(data.X(:)) max(data.X(:))];
            yLim =[min(min(min(data.bounds(1,:,1:nBands,1:nBands2)))) ...
                max(max(max(data.bounds(2,:,1:nBands,1:nBands2))))];
            set(gca, 'LineWidth', 1.5, obj.sfont{:});
            if diff(xLim)>0 && diff(yLim)>0
                set(gca,'XLim',xLim,'yLim',yLim); 
            else
                axis auto
            end
                
            % Add arrows for cross-correlation graphs
            if min(data.X(:))<0
                pos = get(gca,'Position');
                annotation('arrow',[pos(1)+pos(3)/2-pos(3)/100 pos(1)+pos(3)/100],...
                    [pos(2)+pos(4)/100 pos(2)+pos(4)/100],'Linewidth',2);
                annotation('textbox',[pos(1)+pos(3)/10 pos(2)+pos(4)/100 ...
                    pos(3)/2 pos(4)/20],'String',['After ' obj.Input1],'EdgeColor','none',obj.sfont{:});
                annotation('arrow',[pos(1)+pos(3)/2+pos(3)/100 pos(1)+pos(3)-pos(3)/100],...
                    [pos(2)+pos(4)/100 pos(2)+pos(4)/100],'Linewidth',2);
                annotation('textbox',[pos(1)+6*pos(3)/10 pos(2)+pos(4)/100 ...
                    pos(3)/2 pos(4)/20],'String',['Before ' obj.Input1],'EdgeColor','none',obj.sfont{:});
            end
            
            % Create checkboxes if multiple bands
            axesPos = get(get(h(1),'Parent'),'Position');
            mainFig = get(get(h(1),'Parent'),'Parent');
            obj.bands=[];
            if nBands>1 || nBands2>1
                set(mainFig,'Toolbar','figure');
                axesRightPos = axesPos(1)+axesPos(3);
                bandPanel = uipanel(mainFig,...
                    'BackgroundColor',get(mainFig,'Color'),....
                    'BorderType','none',....
                    'Position',[axesRightPos axesPos(2) 1-axesRightPos axesPos(4)]);
                for j=1:nBands2
                    for i=1:nBands
                        ind=sub2ind([nBands nBands2],i,j);
                        obj.bands(ind) = uicontrol(bandPanel,'Style','checkbox',...
                            'BackgroundColor',get(mainFig,'Color'),....
                            'Units','normalized','String',['(' num2str(i) ',' num2str(j) ')'],...
                            'ForegroundColor',colors(ind,:),'Value',1,obj.sfont{:},...
                            'Position',[0 1-ind/nBandsTot 1 1/nBandsTot],...
                            'Callback',@(hObject,event) updateDraw(obj,h,data));
                    end
                end
            end
        end
        
        function updateDraw(obj,h,data)
            nBands = min(size(data.bounds,2),obj.nBandMax);
            nBands2 = min(size(data.bounds,3),obj.nBandMax);
            if nBands>1 || nBands2>1
                states=logical(arrayfun(@(x) get(x,'Value'),obj.bands));
                set(h(states,:),'Visible','on');
                set(h(~states,:),'Visible','off');
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
            params(4).name='XLabel';
            params(4).validator=@ischar;
            params(5).name='YLabel';
            params(5).validator=@ischar;
            params(6).name='Input1';
            params(6).validator=@ischar;
            params(7).name='sfont';
            params(7).validator=@iscell;
            params(8).name='lfont';
            params(8).validator=@iscell;
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end
end