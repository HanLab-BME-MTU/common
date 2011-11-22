classdef WindowsDisplay < MovieDataDisplay
    %Conrete class for displaying flow
    properties
        Color='r';  
        FaceAlpha=.2;
        showNum=5;
    end
    methods        
        function obj=WindowsDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end

        function h=initDraw(obj,data,tag,varargin)

            h=plotWindows(data,{obj.Color,'FaceAlpha',obj.FaceAlpha},obj.showNum);
            set(h,'Tag',tag); 
        end

        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            delete(h);
            h=plotWindows(data,{obj.Color,'FaceAlpha',obj.FaceAlpha},obj.showNum);
            set(h,'Tag',tag);

        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@ischar;
            params(2).name='FaceAlpha';
            params(2).validator=@isscalar;
            params(3).name='showNum';
            params(3).validator=@isscalar;
        end
        function f=getDataValidator()
            f=@iscell;
        end
    end    
end
