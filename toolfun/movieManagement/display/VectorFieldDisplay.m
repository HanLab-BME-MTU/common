classdef VectorFieldDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='k';  
        vectorScale=1;
        vectorCMap=[];
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
  
            if ~isempty(obj.vectorCMap) 
                nColors = size(obj.vectorCMap,1);
                intensity = sqrt((data(:,3)).^2+(data(:,4)).^2);
                maxIntensity=max(intensity(:));
                vColor = floor(intensity/maxIntensity*(nColors-1)+0.5)+1;
 
                h=arrayfun(@(x) quiver(data(vColor==x, 1),data(vColor==x, 2),...
                    obj.vectorScale*(data(vColor==x,3)),...
                    obj.vectorScale*(data(vColor==x,4)),0,...
                    '-','Autoscale',autoscale,...
                    'Color',obj.vectorCMap(x,:),varargin{:}), 1:nColors);   
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
        end
        function updateDraw(obj,h,data)
            if ~isempty(obj.vectorCMap)
                nColors = size(obj.Color,1);
                intensity = sqrt((data(:,3)).^2+(data(:,4)).^2);
                maxIntensity=max(intensity(:));
                vColor = floor(intensity/maxIntensity*(nColors-1)+0.5)+1;
 
                arrayfun(@(x) set(h(x),'XData',data(vColor==x, 1),...
                    'YData',data(vColor==x, 2),...
                    'UData',obj.vectorScale*data(vColor==x, 3),...
                    'VData',obj.vectorScale*data(vColor==x, 4),...
                    'Color',obj.vectorCMap(x,:)),1:nColors);
            else
                set(h,'XData',data(:,1),'YData',data(:,2),...
                    'UData',obj.vectorScale*data(:, 3),...
                    'VData',obj.vectorScale*data(:, 4));
             end
            
            
        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@(x)ischar(x) ||isvector(x));
            ip.addParamValue('vectorScale',obj.vectorScale,@isscalar);      
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@isnumeric;
        end
    end    
end