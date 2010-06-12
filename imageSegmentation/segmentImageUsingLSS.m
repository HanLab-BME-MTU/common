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

% STEP 1: resize image, phi and domain to a dyadic size

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

% STEP 2: define B3-spline sub-sampled filter for step 8 (using in eq. 13)
x = -2^(h+1)+1:2^(h+1)-1;
    
b3h = (1/6) * x.^3 .* (x >= 0) - (4/6) * (x-1).^3 .* (x-1 >= 0) + ...
    (x - 2).^3 .* (x-2 >= 0) - (4/6) * (x-3).^3 .* (x-3 >= 0) + ...
    (1/6) * (x-4).^3 .* (x-4 >= 0);

b3h = b3h / sum(b3h);

% STEP 3: down sample level set using a iterative bicubic method.
scaledPhi = phi;

for s = 1:h
    scaledPhi = imresize(scaledPhi, .5);
end

% STEP 4: initialize B-spline coefficients from scaled phi
b3SplineCoeffs = b3spline2D(scaledPhi);

% STEP 5: normalize B-spline coefficients (eq. 15)
b3SplineCoeffs = b3SplineCoeffs / max(abs(b3SplineCoeffs(:)));

% STEP 6: compute phi value from B-spline coefficients
phi = ib3spline2D(b3SplineCoeffs,size(phi));

%% ------ Main loop ------- %%

epsilon = 1;

iter = 0;
prevEnergy = 0;
energy = +Inf;

while abs(energy - prevEnergy) > eps && iter < maxIter
    
    % STEP 7: compute image feature (eq. 19)
    
    % Define heavy side and dirac functions
    heavySide = .5 + (1 / pi) * atan(phi / epsilon);
    diracFunc = (1 / (pi * epsilon)) / (1 + (phi / epsilon).^2);
    
    % Compute image average in and out levelset
    muIn = sum(sum(ima .* heavySide)) / sum(heavySide(:));
    muOut = sum(sum(ima .* (1 - heavySide))) / sum(1 - heavySide(:));
   
    dataIn = (ima - muIn).^2;
    dataOut = (ima - muOut).^2;
    
    div = dX / normDPhi + dY / normDPhi;
    w = (dataIn.^2 - dataOut.^2 - nu * div) .* diracFunc;
    w = w / max(abs(w(:)));
    
    % STEP 8: compute energy gradient (eq. 13)
    dJ = conv2(b3h, b3h, w, 'same');
    dJ = dJ(1:2^h:end, 1:2^h:end);
        
    % STEP 9: gradient descent feedback adjustement
    newEnergy = +Inf;
    lambda = 1;
    iter2 = 0;
    
    while energy <= newEnergy && iter2 < 5

        % Compute new B3-spline coefficients
        newB3SplineCoeffs = b3SplineCoeffs - lambda * dJ;
        newB3SplineCoeffs = newB3SplineCoeffs / max(abs(newB3SplineCoeffs(:)));
        
        % Compute new level set
        newPhi = ib3spline2D(newB3SplineCoeffs, size(phi));

        % define heavy side and dirac functions on newPhi
        newHeavySide = .5 + (1 / pi) * atan(newPhi / epsilon);
        newDiracFunc = (1 / (pi * epsilon)) / (1 + (newPhi / epsilon).^2);
        
        % Compute new energy function using newPhi
        dX = b3spline1D(newPhi);
        dY = b3spline1D(newPhi')';
        [dX,dY] = gradient(dX,dY);
        newNormDPhi = sqrt(dX.^2 + dY.^2);        

        newJ = dataIn .* newHeavySide + dataOut .* (1 - newHeavySide) + ...
            (nu / newSize^2) * newNormDPhi .* newDiracFunc;
        
        newEnergy = sum(newJ(:));
        
        lambda = lambda / 1.5;
        
        iter2 = iter2 + 1;
    end
    
    % if iter2 == 5, the main loop will stop
    if iter2 ~= 5
        phi = newPhi;
        b3SplineCoeffs = newB3SplineCoeffs;
    end
    
    % update energies
    prevEnergy = energy;
    energy = newEnergy;
    
    iter = iter + 1;
end

BW = phi >= 0;
