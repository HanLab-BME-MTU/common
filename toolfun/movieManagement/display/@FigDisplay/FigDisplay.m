classdef FigDisplay < MovieDataDisplay
    %Concreate class to display general figure plot
    % Andrew R. Jamieson Mar 2017
    
    properties
        plotFunc = @plot; 
        plotFunParams = {};
    end

    methods
        function obj=FigDisplay(varargin)
            obj@MovieDataDisplay(varargin{:})
        end
        
        function h = initDraw(obj, data, tag, varargin)

            ip =inputParser;
            ip.addParameter('Parent', [], @ishandle);
            ip.parse(varargin{:})
            parent_h = ip.Results.Parent;
            if isempty(obj.plotFunParams)
                h = obj.plotFunc(data);    
            else
                h = obj.plotFunc(data, obj.plotFunParams{:});
            end
            set(h,'Tag', tag);
        end
        function updateDraw(obj, h, data)   
        end  
    end 

    methods (Static)
        function params=getParamValidators()
            params(1).name = 'plotFunc';
            params(1).validator = @(A)validateattributes(A,{'function_handle'},{'nonempty'});
            params(2).name = 'plotFunParams';
            params(2).validator = @iscell;
        end
        function f=getDataValidator()
            f=@isstruct;% (A)validateattributes(A,{'struct'},{'nonempty'});
        end
    end    
end