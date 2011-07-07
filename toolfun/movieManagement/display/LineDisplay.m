classdef LineDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
    properties
        Marker = 'none';
        LineStyle = '-'
        Color='r';        
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
        end
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
        end
        function updateDraw(obj,h,data)
            set(h,'XData',data(:,2),'YData',data(:,1))
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@ischar);
            ip.addParamValue('Marker',obj.Marker,@ischar);
            ip.addParamValue('LineStyle',obj.LineStyle,@ischar);
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
        function status=isOverlay()
            status=true;
        end
    end    
end