classdef FigFileDisplay < MovieDataDisplay
    %Concreate class to display external fig-file
    
    methods
        function obj=FigFileDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        
        function h=initDraw(obj,data,tag,varargin)
            % Plot the image and associate the tag
            s=load(data,'-mat');
            fields=fieldnames(s);
            struct2handle(s.(fields{1}).children,gcf);
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            s=load(data,'-mat');
            fields=fieldnames(s);
            struct2handle(s.(fields(1)).children,h);
        end
        
        function additionalInputParsing(obj,ip)
        end
        
        function setProperties(obj,ip)
            
        end
        
    end 
   
 
    methods (Static)
        function f=dataCheck()
            f=@ischar;
        end
    end    
end