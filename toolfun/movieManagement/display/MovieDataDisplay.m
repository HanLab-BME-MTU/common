classdef MovieDataDisplay < handle
    % Abstract class for displaying MovieData components output
    % Delegates drawing methods to the concrete classes
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
            obj.setProperties(ip);
            
            % Retrieve the axes handle and call the create figure method 
            hAxes = ip.Results.hAxes;
            set(hAxes,'NextPlot','add');
            
            % Get the component handle and call the adapted draw function
            h = findobj(hAxes,'-regexp','Tag',tag);
            if ~isempty(h) && ishandle(h)
                obj.updateDraw(h,data);
            else
                h=obj.initDraw(data,tag,'Parent',hAxes);
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