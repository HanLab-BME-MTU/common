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
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
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
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.vectorScale=ip.Results.vectorScale;
            obj.Colormap=ip.Results.Colormap;
%             obj.CLim= ip.Results.CLim;
        end
        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            delete(h);
            obj.initDraw(data,tag);
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@(x)ischar(x) ||isvector(x));
            ip.addParamValue('vectorScale',obj.vectorScale,@isscalar); 
            ip.addParamValue('Colormap',obj.Colormap,@(x) isempty(x) || isnumeric(x));
%             ip.addParamValue('CLim',obj.CLim,@(x) isempty(x) || isnumeric(x));
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end