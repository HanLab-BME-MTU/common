classdef MovieDataDisplay < handle
    % Abstract class for displaying image processing output
    % Based on the template method
    % Sebastien Besson, July 2011
    
    methods
        function h=draw(obj,data,tag,varargin)
            % Template method to draw a movie data component
            
            % Check input
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'MovieDataDisplay'));
            ip.addRequired('data',obj.dataCheck());
            ip.addRequired('tag',@ischar);
            ip.addParamValue('hAxes',gca,@ishandle);
            obj.additionalInputParsing(ip);
            ip.parse(obj,data,tag,varargin{:});
            set(ip.Results.hAxes,'NextPlot','add');
            
            obj.setProperties(ip);
            % Get the component handle
            h = findobj(ip.Results.hAxes,'-regexp','Tag',tag);

            % Call adapted draw function
            if ~isempty(h) && ishandle(h)
                obj.updateDraw(h,data);
            else
                h=obj.initDraw(data,tag,'Parent',ip.Results.hAxes);
            end
        end
    end
    methods(Abstract)
        setProperties(obj,ip)
        additionalInputParsing(obj,ip)
        initDraw(obj,data,tag,varargin)
        updateDraw(obj,h,data,varargin)
    end
    methods (Static,Abstract)
        dataCheck()
    end           
end