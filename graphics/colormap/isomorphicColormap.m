function colormap = isomorphicColormap(color, cLength, test)
%ISOMORPHICCOLORMAP generates isomorphic colormaps
%
% SYNOPSIS: colormap = isomorphicColormap(color,cLength, test)
%
% INPUT color (opt): 'red', 'green', 'blue', 'b/y' for r, g, b, and
%                    blue/yellow colormaps. Default: 'green'
%                    'bw','gw','rw' for blue, green, red going from white
%                    to color
%       cLength (opt): Number of entries of the colormap. Default: 64
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
def_cLength = 64;

% test input
if nargin < 1 || isempty(color)
    color = def_color;
else
    % no need to check here - we will do that with the switch command below
end
if nargin < 2 || isempty(cLength)
    cLength = def_cLength;
end
if nargin < 3 || isempty(test)
    test = def_test;
end


% select colormap
switch color
    case {'red','r'}
        % change color, saturation a tiny bit, luminosity lots
        cmap = hsl2rgb([[linspace(0.95,1,cLength/2)';linspace(0,0.05,cLength/2)'],...
            linspace(0.7,0.8,cLength)',linspace(0.1,1,cLength)']);
    case {'rw'}
        % same as 'r', but going from white to color
        cmap = hsl2rgb([[linspace(0.95,1,cLength/2)';linspace(0,0.05,cLength/2)'],...
            repmat(1,cLength,1),linspace(1,0.5,cLength)']);
    case {'green','g'}
        % change color, saturation a tiny bit, luminosity lots
        cmap = hsl2rgb([linspace(0.3,0.4,cLength)',linspace(0.7,0.8,cLength)',...
            linspace(0.1,1,cLength)']);
        case {'gw'}
        % same as 'g', but going from white to color
        cmap = hsl2rgb([linspace(0.3,0.4,cLength)',...
            repmat(1,cLength,1),linspace(1,0.5,cLength)']);
    case {'blue','b'}
        % change color, saturation a tiny bit, luminosity lots
        cmap = hsl2rgb([linspace(0.5,0.6,cLength)',linspace(0.7,0.8,cLength)',...
            linspace(0.1,1,cLength)']);
        case {'bw'}
        % same as 'g', but going from white to color
        cmap = hsl2rgb([linspace(0.5,0.6,cLength)',...
            repmat(1,cLength,1),linspace(1,0.5,cLength)']);
    case 'b/y'
        % color: blue->yellow, saturation 1->0->1, luminosity 0.5 (=max
        % color)
        cmap = hsl2rgb([repeatEntries([0.66;0.16],cLength/2),[linspace(0.8,0,cLength/2)';...
        linspace(0,0.8,cLength/2)'],ones(cLength,1)*0.50]);
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