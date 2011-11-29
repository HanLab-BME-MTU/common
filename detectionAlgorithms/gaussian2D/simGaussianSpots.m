function [frame, xv, yv, sv, Av] = simGaussianSpots(nx, ny, sigma, varargin)
% SIMGAUSSIANSPOTS generates a given number of 2D Gaussians in an image.
% The generated Gaussian signals do not overlap with the image boundaries.
%
%   Usage: [frame,xv,yv,sv,Av]=simGaussianSpots(nx,ny,sigma,varargin)
%   Input:
%           nx -> image size in x direction
%           ny -> image size in y direction
%           sigma -> 2D Gaussian standard deviation (in pixels)
%       varargin may include:
%           x -> x coordinates of centers of 2D Gaussians (can be subpixel)
%           y -> y coordinates of centers of 2D Gaussians (can be subpixel)
%           A -> amplitudes of 2D Gaussians
%           npoints -> number of 2D Gaussian to be generated
%           background -> value of background
%
%   Output:
%           frame -> image with 2D Gaussian signals
%           xv -> x coordinates of centers of Gaussians (can be subpixel)
%           vy -> y coordinates of centers of Gaussians (can be subpixel)
%           sv -> vector of standard deviation
%           Av -> vector of amplitudes

% Francois Aguet, last modified June 30 2011

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('nx',@isnumeric);
ip.addRequired('ny',@isnumeric);
ip.addRequired('sigma',@isnumeric);
ip.addParamValue('x', []);
ip.addParamValue('y', []);
ip.addParamValue('A', []);
ip.addParamValue('npoints', 1);
ip.addParamValue('background', 0);
ip.addParamValue('verbose', 'off', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.parse(nx, ny, sigma, varargin{:});

np = ip.Results.npoints;
c = ip.Results.background;
xv = ip.Results.x;
yv = ip.Results.y;
Av = ip.Results.A;

if length(xv) ~= length(yv)
    error('''x'' and ''y'' must be of same size');
end

if ~isempty(xv)
    np=length(xv);
end

% feature vectors
sv = sigma*ones(1,np);
wv = ceil(4*sv);

% discard signals too near to image boundaries
if ~isempty(xv)
    idx=xv > max(wv)+1 & xv < nx-max(wv)-1;
    idy=yv > max(wv)+1 & yv < ny-max(wv)-1;
    id=idx & idy;    
    if strcmpi(ip.Results.verbose, 'on')
        fprintf('Number of discarded points: %d\n', length(xv) - sum(id));
    end
    xv=xv(id);
    yv=yv(id);
    sv=sv(id);
    np=length(xv);
end

w_max = 2*max(wv)+1;
if (w_max>nx || w_max>ny)
    error('Dimensions are too small.');
end

if isempty(xv)
    xv = (nx-2*wv-1).*rand(1,np) + wv+1;
    yv = (ny-2*wv-1).*rand(1,np) + wv+1;
elseif length(xv)~=length(yv)
    error('''x'' and ''y'' must be of same length.');
end

if isempty(Av)
    Av = ones(1,np);
end

frame = c*ones(ny, nx);

[~, idx] = sort(xv+yv*nx); % sort spots according to row position
xv = xv(idx);
yv = yv(idx);
Av = Av(idx);

for k = 1:np
    xi = round(xv(k));
    yi = round(yv(k));
    dx = xv(k)-xi;
    dy = yv(k)-yi;
    wi = wv(k);
    [x,y] = meshgrid(-wi:wi);
    r2 = (x-dx).^2+(y-dy).^2;
    g = Av(k) * exp(-r2/(2*sv(k)^2));
    frame(yi-wi:yi+wi,xi-wi:xi+wi) = frame(yi-wi:yi+wi,xi-wi:xi+wi) + g;
end
