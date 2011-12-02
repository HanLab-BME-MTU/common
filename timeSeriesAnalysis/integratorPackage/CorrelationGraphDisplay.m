classdef CorrelationGraphDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        LineStyle = '-';
        LineWidth = 2;
        Color='r';
        Input1 ='';
        Input2 ='';
    end
    methods
        function obj=CorrelationGraphDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
        function h=initDraw(obj,data,tag,varargin)
            
            % define small and large fonts
            sfont = {'FontName', 'Helvetica', 'FontSize', 18};
            lfont = {'FontName', 'Helvetica', 'FontSize', 22};
            
            [nx,nBands] = size(data.avgCorrFun);
            h=-ones(nBands,2);
            colors = hsv(nBands);
            for i=1:nBands
                h(i,1)=plot(data.lags(:,i),data.avgCorrFun(:,i),'Line',obj.LineStyle,...
                    'LineWidth',obj.LineWidth,'Color',colors(i,:));
                hold on
                h(i,2)=errorbar(data.lags(:,i),data.avgCorrFun(:,i),data.steCorrFun(:,i),...
                    'LineWidth', 2,'Color',colors(i,:));
            end
            xLim=[min(data.lags(:)) max(data.lags(:))];
            yLim =[min(data.avgCorrFun(:)-data.steCorrFun(:)) max(data.avgCorrFun(:)+data.steCorrFun(:))];
            xlabel('Lag (s)',lfont{:})
            ylabel('Correlation',lfont{:})
            set(gca, 'LineWidth', 1.5, sfont{:},'XLim',xLim,'YLim',yLim);
            
            if ~isempty(data.avgBounds)
                
                upline  = repmat(data.avgBounds(1,:),nx,1);
                h(3)=plot(data.lags,upline,'Linewidth',2,'Color','r');
                
                dline  = repmat(data.avgBounds(2,:),nx,1);
                h(4)=plot(data.lags,dline,'Linewidth',2,'Color','r');
            end
            
            if min(data.lags)<0
                pos = get(gca,'Position');
                annotation('arrow',[pos(1)+pos(3)/2-pos(3)/100 pos(1)+pos(3)/100],...
                    [pos(2)+pos(4)/100 pos(2)+pos(4)/100],'Linewidth',2);
                annotation('textbox',[pos(1)+pos(3)/10 pos(2)+pos(4)/100 ...
                    pos(3)/2 pos(4)/20],'String',['After ' obj.Input1],'EdgeColor','none',sfont{:});
                annotation('arrow',[pos(1)+pos(3)/2+pos(3)/100 pos(1)+pos(3)-pos(3)/100],...
                    [pos(2)+pos(4)/100 pos(2)+pos(4)/100],'Linewidth',2);
                annotation('textbox',[pos(1)+6*pos(3)/10 pos(2)+pos(4)/100 ...
                    pos(3)/2 pos(4)/20],'String',['Before ' obj.Input1],'EdgeColor','none',sfont{:});
            end
            
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            axes(get(get(h(1),'Parent')));
            tag = get(h(1),'Tag');
            cla
            obj.initDraw(data,tag);
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
            params(4).name='Input1';
            params(4).validator=@ischar;
            params(5).name='Input2';
            params(5).validator=@ischar;
        end
        
        function f=getDataValidator()
            f=@isstruct;
        end
    end
end