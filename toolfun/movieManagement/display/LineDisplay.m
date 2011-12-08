classdef LineDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Color='r';
        Marker = 'none';
        LineStyle = '-';
        LineWidth = 1;
        XLabel='';
        YLabel='';
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
    end
    methods
                
        function obj=LineDisplay(varargin)
            obj@MovieDataDisplay(varargin{:})
        end
        function h=initDraw(obj,data,tag,varargin)
            
            h=plot(data(:,2),data(:,1),varargin{:});
            set(h,'Tag',tag,'Color',obj.Color,'Marker',obj.Marker,...
                'Linestyle',obj.LineStyle,'LineWidth',obj.LineWidth);
            obj.setAxesProperties;
            
        end
        function updateDraw(obj,h,data)
            set(h,'XData',data(:,2),'YData',data(:,1));
            obj.setAxesProperties;
        end
    end
    methods(Access=protected)
        function setAxesProperties(obj)
            if ~isempty(obj.XLabel),xlabel(obj.XLabel,obj.lfont{:}); end
            if ~isempty(obj.YLabel),ylabel(obj.YLabel,obj.lfont{:}); end
            set(gca,'LineWidth', 1.5, obj.sfont{:})
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
            params(4).name='LineWidth';
            params(4).validator=@isscalar;
            params(5).name='XLabel';
            params(5).validator=@ischar;
            params(6).name='YLabel';
            params(6).validator=@ischar;
            params(7).name='sfont';
            params(7).validator=@iscell;
            params(8).name='lfont';
            params(8).validator=@iscell;
        end
        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end