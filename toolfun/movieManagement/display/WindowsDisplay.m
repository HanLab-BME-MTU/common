classdef WindowsDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='r';  
        FaceAlpha=.2;
        showNum=5;
    end
    methods
        function obj=WindowsDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)

            h=plotWindows(data,{obj.Color,'FaceAlpha',obj.FaceAlpha},obj.showNum);
            set(h,'Tag',tag); 
        end
        function setProperties(obj,ip)
            obj.Color=ip.Results.Color;
            obj.FaceAlpha=ip.Results.FaceAlpha;
        end
        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            delete(h);
            h=plotWindows(data,{obj.Color,'FaceAlpha',obj.FaceAlpha},obj.showNum);
            set(h,'Tag',tag);

        end
        function additionalInputParsing(obj,ip)
            ip.addParamValue('Color',obj.Color,@(x)ischar(x) ||isvector(x));
            ip.addParamValue('FaceAlpha',obj.FaceAlpha,@isscalar);
            ip.addParamValue('showNum',obj.showNum,@isscalar);
        end 
    end    
    
    methods (Static)
        function f=dataCheck()
            f=@iscell;
        end
    end    
end
