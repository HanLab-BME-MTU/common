classdef VectorFieldDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='r';        
    end
    methods
        function h=initDraw(obj,data,tag,varargin)
            h=quiver(data(:, 2),data(:, 1),data(:, 4)-data(:, 2),...
                data(:, 3)-data(:, 1),'-','Color',obj.Color,varargin{:});
            set(h,'Tag',tag);
        end
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
        end
        function updateDraw(obj,h,data)
            set(h,'XData',data(:,2),'YData',data(:,1),...
                'UData',data(:, 4)-data(:, 2),'VData',data(:, 3)-data(:, 1))
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@ischar);
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