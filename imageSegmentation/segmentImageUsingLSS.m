function [BW,phi,iter] = segmentImageUsingLSS(ima, domain, phi, nu, h, maxIter)

%% ------ Parse inputs --------- %%

if nargin ~= 6
    error('Invalid number of input argument.');
end

if h ~= floor(h) || h < 0
    error('scale parameter h needs to be a positive integer.');
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
    radius = (1/3) * min(nrows,ncols);
    [X Y] = meshgrid((1:ncols) - cx, (1:nrows) - cy);
    phi = -sqrt(X.^2 + Y.^2) + radius;
    
elseif size(phi,1) ~= nrows || size(phi,2) ~= ncols
    error('ima size and phi size differ.');
end

%% ------ Initialization ------ %%

% STEP 1: resize image, phi and domain to a squared power of 2 size width

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

% STEP 2: down sample level set using a iterative bicubic method.

scaledPhi = phi;

for s = 1:h
    scaledPhi = imresize(scaledPhi, .5);
end

% STEP 3: initialize B-spline coefficients from scaled phi

BS = convertRowSignalToBSplineCoefficients(scaledPhi);
BS = convertRowSignalToBSplineCoefficients(BS');
BS = BS';

% STEP 4: normalize B-spline coefficients according to eq. 15

BS = BS / max(abs(BS(:)));

% STEP 5: update level set value from down-sampled B-spline coefficients
% using B3-spline interpolation

[X,Y] = meshgrid(1:newSize,1:newSize);

% Rescale levelset coordinates into the B-spline space (side with = 2^(p-h))
xScaled = X / 2^h;
yScaled = Y / 2^h;

% The subsequent code can be replaced by the following line (TO BE CHECKED)
%phi = interp2(X / 2^h, Y / 2^h, bSplineCoeffs, X, Y, 'spline');

% Compute the interpolation indexes
xIndex = arrayfun(@(i) floor(xScaled) + i, -1:2, 'UniformOutput', false);
xIndex = cat(3,xIndex{:});

yIndex = arrayfun(@(i) floor(yScaled) + i, -1:2, 'UniformOutput', false);
yIndex = cat(3,yIndex{:});

% Compute weights
xW = xScaled - xIndex(:,:,2);
yW = yScaled - yIndex(:,:,2);

xWeight = zeros([size(xScaled), 4]);
xWeight(:,:,4) = (1/6) * xW.^3;
xWeight(:,:,1) = (1/6) + (1/2) * xW .* (xW - 1) - xWeight(:,:,4);
xWeight(:,:,3) = xW + xWeight(:,:,1) - 2 * xWeight(:,:,4);
xWeight(:,:,2) = 1 - xWeight(:,:,1) - xWeight(:,:,3) - xWeight(:,:,4);

yWeight = zeros([size(yScaled), 4]);
yWeight(:,:,4) = (1/6) * yW.^3;
yWeight(:,:,1) = (1/6) + (1/2) * yW .* (yW - 1) - yWeight(:,:,4);
yWeight(:,:,3) = yW + yWeight(:,:,1) - 2 * yWeight(:,:,4);
yWeight(:,:,2) = 1 - yWeight(:,:,1) - yWeight(:,:,3) - yWeight(:,:,4);

% Add mirror condition at border 
BSpadded = padarrayXT(BS, [2 2], 'symmetric', 'both');

% Compute interpolation for each X,Y = 
% sum_k=1^4 yWeight(k) sum_l=1^4 xWeight(l) BSpadded(yIndex(l), xIndex(l))

phi(:) = 0;
tmp = zeros(2^p);

for k = 1:4
    tmp(:) = 0;
    
    yShifted = yIndex(:,:,k) + 2;
    
    for l = 1:4
        xShifted = xIndex(:,:,l) + 2;
        
        ind = sub2ind(size(BSpadded), yShifted(:), xShifted(:));
        
        tmp = tmp + xWeight(:,:,l) .* reshape(BSpadded(ind), [512, 512]);
        
    end
    
    phi = phi + yWeight(:,:,k) .* tmp;
end


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

for i = 2:n
    in(:,i) = in(:,i) + z * in(:,i-1);
end

in(:,end) = (z * in(:,end-1) + in(:,end)) * z / (z * z - 1);

for i = n-1:-1:1
    in(:,i) = z * (in(:,i+1) - in(:,i));
end
