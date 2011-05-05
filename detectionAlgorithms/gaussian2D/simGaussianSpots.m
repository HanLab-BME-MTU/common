% Francois Aguet, last modified April 20 2011

function [frame, xv, yv, sv, Av] = simGaussianSpots(nx, ny, sigma, varargin)

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
ip.parse(nx, ny, sigma, varargin{:});

np = ip.Results.npoints;
c = ip.Results.background;
xv = ip.Results.x;
yv = ip.Results.y;
Av = ip.Results.A;

% feature vectors
sv = sigma*ones(1,np);
wv = ceil(4*sv);

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
