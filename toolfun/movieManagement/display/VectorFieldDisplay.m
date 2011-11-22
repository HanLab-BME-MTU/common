classdef VectorFieldDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='k';  
        vectorScale=1;
        Colormap=[];
        CLim = [];
    end
    methods
        function obj=VectorFieldDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
        
        function h=initDraw(obj,data,tag,varargin)
            if isempty(obj.vectorScale),autoscale='on'; else autoscale='off'; end
  
            if ~isempty(obj.Colormap) 
                nColors = size(obj.Colormap,1);
                intensity= (data(:,3).^2+data(:, 4).^2).^(1/2);
                if ~isempty(obj.CLim)
                    vColor = floor((intensity-obj.CLim(1))/diff(obj.CLim)*(nColors-1)+0.5)+1;
                else
                    vColor = floor(intensity/max(intensity(:))*(nColors-1)+0.5)+1;
                end
 
                h=arrayfun(@(x) quiver(data(vColor==x, 1),data(vColor==x, 2),...
                    obj.vectorScale*(data(vColor==x,3)),...
                    obj.vectorScale*(data(vColor==x,4)),0,...
                    '-','Autoscale',autoscale,...
                    'Color',obj.Colormap(x,:),varargin{:}), unique(vColor));   
            else
                
                h=quiver(data(:,1),data(:, 2),obj.vectorScale*(data(:,3)),...
                    obj.vectorScale*(data(:,4)),0,'-','Autoscale',autoscale,...
                    'Color',obj.Color,varargin{:});
         
            end
               
            set(h,'Tag',tag);
        end

        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            delete(h);
            obj.initDraw(data,tag);
        end
    end    
    
    methods (Static)
         function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@(x)ischar(x) ||isvector(x);
            params(2).name='vectorScale';
            params(2).validator=@isscalar;
            params(3).name='Colormap';
            params(3).validator=@(x) ischar(x) || isnumeric(x);
            params(4).name='CLim';
            params(4).validator=@isvector;
        end

        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end