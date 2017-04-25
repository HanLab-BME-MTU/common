function [ hpts ] = linkaxes_pointer( axes_handles, option, varargin )
%linkaxes_pointer Links axes and creates a common impoint for all axes
%
% axes_handles: array of axes handles, default: all axes
% option: 'x', 'y', 'xy', or 'off'. Use 'off' if you just want the pointers
%
% See also linkaxes

% Parse extra parameters
ip = inputParser;
ip.addParameter('syncColors',false,@islogical);
ip.addParameter('drawCrosshairs',false,@islogical);
ip.parse(varargin{:});

if(nargin < 1 || isempty(axes_handles))
    axes_handles = [];
end
if(nargin < 2 || isempty(option))
    option = 'xy';
end

lines = [];

if(isempty(axes_handles))
    axes_handles = findobj('-class','matlab.graphics.axis.Axes');
end

linkaxes(axes_handles,option);

hpts = arrayfun(@createPointer,axes_handles,'UniformOutput',false);
hpts = [hpts{:}];

if(ip.Results.syncColors)
    linkprop(lines,'Color');
end

    function hpt = createPointer(ax)
        hpt = impoint(ax,mean([xlim(ax); ylim(ax)],2)');
        addNewPositionCallback(hpt,@syncPositions);
        if(ip.Results.drawCrosshairs)
            plus = findobj(hpt,'Tag','plus');
            set(plus,'MarkerSize',100);
        end
        if(ip.Results.syncColors)
            lines = [lines; findobj(hpt,'Type','Line')];
        end
    end
    function syncPositions(pos)
        for pt = 1:length(hpts)
            if(isvalid(hpts(pt)))
                setPosition(hpts(pt),pos);
            end
        end
    end
    function syncColors(h, ev)
    end

end

