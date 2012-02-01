% Francois Aguet, 02/01/2012

function [res, theta, nms, scaleMap] = multiscaleSteerableDetector(img, M, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addRequired('M', @(x) ismember(x, 1:5));
ip.addOptional('sigmaArray', [1 2 4]);
ip.parse(img, M, varargin{:});

sigma = ip.Results.sigmaArray;
ns = numel(sigma);

% arrays to store results from individual scales
ires = cell(1,ns);
itheta = cell(1,ns);
inms = cell(1,ns);

% responses at individual scales
for si = 1:ns
    [ires{si}, itheta{si}, inms{si}] = steerableDetector(img, ip.Results.M, sigma(si));
    ires{si} = sigma(si) * ires{si}; % scale normalization
end

% determine  f(x,y) = argmax_s r_s(x,y) (init. at scale 1)
res = ires{1};
theta = itheta{1};
nms = inms{1};
scaleMap = ones(size(img));

for si = 2:ns
    idx = ires{si}>res;
    res(idx) = ires{si}(idx);
    nms(idx) = inms{si}(idx);
    theta(idx) = itheta{si}(idx);
    scaleMap(idx) = si;
end