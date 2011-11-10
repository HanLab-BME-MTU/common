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
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            h=plot(data(:,2),data(:,1),varargin{:});
            set(h,'Tag',tag,'Color',obj.Color,'Marker',obj.Marker,...
                'Linestyle',obj.LineStyle);
            
            if ~isempty(obj.XLabel),xlabel(obj.XLabel,varargin{:}); end
            if ~isempty(obj.YLabel),ylabel(obj.YLabel,varargin{:}); end
            
        end
        function updateDraw(obj,h,data)
            set(h,'XData',data(:,2),'YData',data(:,1))
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@ischar);
            ip.addParamValue('Marker',obj.Marker,@ischar);
            ip.addParamValue('LineStyle',obj.LineStyle,@ischar);  
            ip.addParamValue('XLabel',obj.XLabel,@ischar); 
            ip.addParamValue('YLabel',obj.YLabel,@ischar); 
        end 
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.Marker=ip.Results.Marker;
            obj.LineStyle=ip.Results.LineStyle;
            obj.XLabel=ip.Results.XLabel;
            obj.YLabel=ip.Results.YLabel;
        end
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end