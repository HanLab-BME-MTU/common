function [BW,phi,iter] = segmentImageUsingLSS(ima, domain, phi, nu, h, maxIter)

%% ------ Parse inputs --------- %%

if nargin ~= 5
    error('Invalid number of input argument.');
end

if ~isinteger(h)
    error('scale parameter h needs to be an integer.');
end

[nrows,ncols] = size(ima);
p = nextpow2(max(nrows,ncols));

if h >= p + 1
    error('scale parameter h greater than image size.');
end

if isempty(domain)
    domain = true(nrows,ncols);
    
elseif size(domain,1) ~= nrows || size(domain,2) ~= ncols
    error('ima size and domain size differ.');
end

if isempty(phi)
    % Create an initial level set from a cone function
    cx = floor(ncols/2);
    cy = floor(nrows/2);
    radius = (2/3) * min(nrows,ncols);
    [X Y] = meshgrid(1:ncols - cx, 1:nrows - cy);
    phi = -sqrt(X.^2 + Y.^2 - radius.^2);
    
elseif size(phi,1) ~= nrows || size(phi,2) ~= ncols
    error('ima size and phi size differ.');
end

%% ------ Initialization ------ %%

% resize image, phi and domain to a squared power of 2 size

newSize = 2^p;

tmp = ima;
ima = zeros(newSize);
ima(1:nrows,1:ncols) = tmp;

tmp = domain;
domain = logical(newSize);
domain(1:nrows,1:ncols) = tmp;
validDomainIdx = find(domain == true);

tmp = phi;
phi = zeros(newSize);
phi(1:nrows,1:ncols) = tmp;

% down sample level set using a iterative bicubic method.

scaledPhi = phi;

for s = 1:h
    scaledPhi = imresize(scaledPhi, .5);
end

% Initialize B-spline coefficients from scaled phi

bSplineCoeffs = convertRowSignalToBSplineCoefficients(scaledPhi);
bSplineCoeffs = convertRowSignalToBSplineCoefficients(bSplineCoeffs');
bSplineCoeffs = bSplineCoeffs';

% Normalize B-spline coefficients according to eq. 15

bSplineCoeffs = bSplineCoeffs / max(abs(bSplineCoeffs(:)));

% Update level set from B-spline coefficients


%% ------ Main loop ------- %%

epsilon = 1;

iter = 0;
prevEnergy = +inf;
energy = 0;

while abs(energy - prevEnergy) > eps && iter < maxIter
    
    prevEnergy = energy;
    
    % define heavy side function
    heavySide = .5 + (1 / pi) * atan(phi / epsilon);
    
    % define dirac function
    diracFunc = (1 / (pi * epsilon)) / (1 + (phi / epsilon).^2);
    
    % Compute image average in and out levelset
    muIn = sum(sum(ima .* heavySide)) / sum(heavySide(:));
    muOut = sum(sum(ima .* (1 - heavySide))) / sum(1 - heavySide(:));
    
    % Compute gradient magnitude of phi
    % TODO
    
    energyFunc = (ima - muIn).^2 .* heavySide + (ima - muOut).^2 .* (1 - heavySide) + ...
        (nu / newSize^2) * normDPhi .* diracFunc;
    
    energy = sum(energyFunc(:));
    
    
    iter = iter + 1;
end


function in = convertRowSignalToBSplineCoefficients(in)

n = size(in,1);

z = sqrt(3) - 2;

in = in * (1 - z) * (1 - 1 / z);

z1 = z;
zn = z^(n-1);
sum = in(:,1) + zn * in(:,end);
horizon = min(n, 2 + floor(log(1e-6) / log(abs(z))));

zn = zn * zn;

for i = 1:horizon-2
    zn = zn / z;
    sum = sum + (z1 + zn) * in(:,i);
    z1 = z1 * z;
end

in(:,1) = sum / (1 - z^(2 * n - 2));

in(:,2:end) = in(:,2:end) + z * in(:,1:end-1);

in(:,end) = (z * in(:,end-1) + in(:,end)) * z / (z * z - 1);

in(:,end-1:-1:1) = z * (in(:,end:-1:2) - in(:,end-1:-1:1));
