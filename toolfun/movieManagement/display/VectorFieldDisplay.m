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
            if isempty(obj.scale),autoscale='on'; else autoscale='off'; end
            h=quiver(data(:, 2),data(:, 1),obj.scale*(data(:, 4)-data(:, 2)),...
                obj.scale*(data(:, 3)-data(:, 1)),0,'-','Autoscale',autoscale,...
                'Color',obj.Color,varargin{:});
            set(h,'Tag',tag);
        end
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.scale=ip.Results.scale;
        end
        function updateDraw(obj,h,data)
            set(h,'XData',data(:,2),'YData',data(:,1),...
                'UData',obj.scale*(data(:, 4)-data(:, 2)),...
                'VData',obj.scale*(data(:, 3)-data(:, 1)))
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@(x)ischar(x) ||isvector(x));
            ip.addParamValue('scale',obj.scale,@isscalar);
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end