classdef VectorFieldDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='k';  
        scale=1;
    end
    methods
        function obj=VectorFieldDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
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
            ip.addParamValue('Color',obj.Color,@(x)ischar(x) ||isvector(x));
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end