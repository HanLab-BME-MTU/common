classdef LineDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        LineStyle = '-'
        Color='r';
        XLabel='';
        YLabel='';
    end
    methods
                
        function obj=LineDisplay(varargin)
            obj@MovieDataDisplay(varargin{:})
        end
        function h=initDraw(obj,data,tag,varargin)
            
            % Set the fonts
            sfont = {'FontName', 'Helvetica', 'FontSize', 18};
            lfont = {'FontName', 'Helvetica', 'FontSize', 22};
            
            h=plot(data(:,2),data(:,1),varargin{:});
            set(h,'Tag',tag,'Color',obj.Color,'Marker',obj.Marker,...
                'Linestyle',obj.LineStyle);
            
            if ~isempty(obj.XLabel),xlabel(obj.XLabel,lfont{:},varargin{:}); end
            if ~isempty(obj.YLabel),ylabel(obj.YLabel,lfont{:},varargin{:}); end
            set(gca,'LineWidth', 1.5, sfont{:})
        end
        function updateDraw(obj,h,data)
            set(h,'XData',data(:,2),'YData',data(:,1))
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
        end
        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end