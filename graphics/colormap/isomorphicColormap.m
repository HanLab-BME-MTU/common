function colormap = isomorphicColormap(color, test)
%ISOMORPHICCOLORMAP generates isomorphic colormaps
%
% SYNOPSIS: colormap = isomorphicColormap(color, test)
%
% INPUT color (opt): 'red', 'green', 'blue', 'b/y' for r, g, b, and
%                    blue/yellow colormaps. Default: 'green'
%       test (opt) : if 1, colormapTest is being run with the colormap.
%                    Default: 0
% OUTPUT colormap: colormap (64x3 array of RGB values)
%			
%
% REMARKS see
%           http://www.research.ibm.com/people/l/lloydt/color/color.htm
%         and
%           http://www.research.ibm.com/dx/proceedings/pravda/index.htm
%         for details
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 02-Apr-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defaults
def_color = 'red';
def_test = false;

% test input
if nargin < 1 || isempty(color)
    color = def_color;
else
    % no need to check here - we will do that with the switch command below
end
if nargin < 2 || isempty(test)
    test = def_test;
end


% select colormap
switch color
    case 'red'
        % change color, saturation a tiny bit, luminosity lots
        cmap = hsl2rgb([[linspace(0.95,1,32)';linspace(0,0.05,32)'],...
            linspace(0.7,0.8,64)',linspace(0.1,1,64)']);
    case 'green'
        % change color, saturation a tiny bit, luminosity lots
        cmap = hsl2rgb([linspace(0.3,0.4,64)',linspace(0.7,0.8,64)',...
            linspace(0.1,1,64)']);
    case 'blue'
        % change color, saturation a tiny bit, luminosity lots
        cmap = hsl2rgb([linspace(0.5,0.6,64)',linspace(0.7,0.8,64)',...
            linspace(0.1,1,64)']);
    case 'b/y'
        % color: blue->yellow, saturation 1->0->1, luminosity 0.5 (=max
        % color)
        cmap = hsl2rgb([repeatEntries([0.66;0.16],32),[linspace(0.8,0,32)';...
        linspace(0,0.8,32)'],ones(64,1)*0.50]);
    otherwise
        error('color %s not implemented yet',color)
end

% check test
if test
    colormapTest(cmap);
end

% check for output argument
if nargout > 0
    colormap = cmap;
end