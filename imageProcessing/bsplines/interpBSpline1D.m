%[fi, fi_dx, fi_d2x] = interpBSpline1D(xi, f, varargin)
%
% Inputs:
%           xi : interpolation coordinates
%            f : function to interpolate
%
% Optional inputs:
%            n : degree of the spline, 1-3. Default: 3
%     boundary : boundary conditions: 'symmetric' (default) or 'periodic'
%
% Outputs:
%           fi : interpolated values at xi
%        fi_dx : derivative of f at xi 
%       fi_d2x : 2nd derivative of f at xi

% Francois Aguet, 2009 (Last modified: 10/05/2011)

function [fi, fi_dx, fi_d2x] = interpBSpline1D(xi, f, varargin)

ip = inputParser;
ip.addRequired('xi');
ip.addRequired('f', @isvector);
ip.addOptional('n', 3, @(x) ismember(x, 1:3));
ip.addOptional('boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic'})));
ip.parse(xi, f, varargin{:});

c = computeBSplineCoefficients(f, ip.Results.n, ip.Results.boundary);

fi = arrayfun(@(i) interpBSplineValue(i, c, ip.Results.n, ip.Results.boundary), xi);
if nargout>1
    fi_dx = arrayfun(@(i) interpBSplineValue(i+0.5, c, ip.Results.n-1, ip.Results.boundary) - ...
        interpBSplineValue(i-0.5, c, ip.Results.n-1, ip.Results.boundary), xi);
end
if nargout>2
    fi_d2x = arrayfun(@(i) interpBSplineValue(i+1, c, ip.Results.n-2, ip.Results.boundary) - ...
        2*interpBSplineValue(i, c, ip.Results.n-2, ip.Results.boundary) + ...
        interpBSplineValue(i-1, c, ip.Results.n-2, ip.Results.boundary), xi);
end
