function [phi,mask,iter] = evolveSplineB3LevelSet(ima, domain, phi, nu, h, maxIter)

%% ------ Parse inputs --------- %%

ima = double(ima);

if nargin ~= 6
    error('Invalid number of input argument.');
end

if h ~= floor(h) || h < 0
    error('scale parameter h needs to be a positive integer.');
end

[nrows,ncols] = size(ima);
p = nextpow2(min(nrows,ncols));

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

% newSize = 2^p;
% 
% tmp = ima;
% ima = zeros(newSize);
% ima(1:nrows,1:ncols) = tmp;
% 
% tmp = domain;
% domain = false(newSize);
% domain(1:nrows,1:ncols) = tmp;
ind = find(domain == true);
% 
% tmp = phi;
% phi = zeros(newSize);
% phi(1:nrows,1:ncols) = tmp;

% STEP 2: define B3-spline sub-sampled filter for step 8 (using in eq. 13)
x = -2:2^-h:2;
x = x(2:end-1);
b3h = sum(cell2mat(arrayfun(@(j) ((-1)^j / 6) * nchoosek(4,j) * ...
    (x + 2 - j).^3 .* (x + 2 - j >= 0), (0:4)', 'UniformOutput', false)));
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

% max iteration for the Gradient Descent Feedback Adjustment (GDFA)
maxIterGDFA = 5;

epsilon = 1;
iter = 0;
prevEnergy = +Inf;
energy = +Inf;

% hFig = figure;
% imagesc(ima); axis image; colormap gray; axis off; hold on;

while energy <= prevEnergy  && energy > eps && iter < maxIter
    
%     % DEBUG
%     hLayers = findall(get(gca, 'Children'), 'Tag', 'contour');
%     delete(hLayers);
%     edges = edge(phi > 0);
%     [y x] = ind2sub(size(phi),find(edges));
%     line(x,y,'LineStyle','none','Marker','.','Tag', 'contour');
%     refresh; figure(hFig);
%     % END OF DEBUG
    
    % STEP 7: compute image feature (eq. 19)
    heavySide = .5 + (1 / pi) * atan(phi / epsilon);
    diracFunc = (1 / (pi * epsilon)) ./ (1 + (phi / epsilon).^2);
    muIn = sum(sum(ima(ind) .* heavySide(ind))) / sum(heavySide(ind));
    muOut = sum(sum(ima(ind) .* (1 - heavySide(ind)))) / sum(1 - heavySide(ind));
    dataIn = (ima - muIn).^2;
    dataOut = (ima - muOut).^2;
    dX = gradient(b3spline1D(phi));
    dY = gradient(b3spline1D(phi'))';
    normDPhi = sqrt(dX.^2 + dY.^2) + eps;
    div = dX ./ normDPhi + dY ./ normDPhi;
    w = (dataIn - dataOut - nu * div) .* diracFunc;
    w = w / max(abs(w(:)));
    
    % STEP 8: compute energy gradient (eq. 13)
    dJ = conv2(b3h, b3h, w, 'same');
    dJ = dJ(1:2^h:end, 1:2^h:end);
        
    % STEP 9: gradient descent feedback adjustement
    newEnergy = +Inf;
    lambda = 1;
    iter2 = 0;
    
    while energy <= newEnergy && iter2 < maxIterGDFA

        % Compute new B3-spline coefficients
        newB3SplineCoeffs = b3SplineCoeffs - lambda * dJ;
        newB3SplineCoeffs = newB3SplineCoeffs / max(abs(newB3SplineCoeffs(:)));
        
        % Compute new level set
        newPhi = ib3spline2D(newB3SplineCoeffs, size(phi));

        % define heavy side and dirac functions on newPhi
        newHeavySide = .5 + (1 / pi) * atan(newPhi / epsilon);
        newDiracFunc = (1 / (pi * epsilon)) ./ (1 + (newPhi / epsilon).^2);
        
        % Compute new energy function using newPhi
        dX = gradient(b3spline1D(newPhi));
        dY = gradient(b3spline1D(newPhi'))';
        newNormDPhi = sqrt(dX.^2 + dY.^2);

        newJ = dataIn .* newHeavySide + dataOut .* (1 - newHeavySide) + ...
            nu * newNormDPhi .* newDiracFunc;
        
        newEnergy = sum(newJ(:));
        
        lambda = lambda / 1.5;
        
        iter2 = iter2 + 1;
    end
    
    % if iter2 == maxIterGDFA, the main loop will stop
    if iter2 ~= maxIterGDFA
        phi = newPhi;
        b3SplineCoeffs = newB3SplineCoeffs;
    end
    
    % update energies
    prevEnergy = energy;
    energy = newEnergy;
    
    iter = iter + 1;
end

% close(hFig);

mask = phi > 0;

% add border
mask = padarray(mask, [1 1], 'replicate');
ind = find(mask(1,:));
if ~isempty(ind)
    mask(1, 1:ind(end)) = true;
end

ind = find(mask(end,:));
if ~isempty(ind)
    mask(end, 1:ind(end)) = true;
end

ind = find(mask(:, 1));
if ~isempty(ind)
    mask(1:ind(end), 1) = true;
end

ind = find(mask(:, end));
if ~isempty(ind)
    mask(1:ind(end), end) = true;
end

mask = imfill(mask,'holes');

mask = mask(2:end-1,2:end-1);

% select the biggest area
s = regionprops(mask,'Area','PixelIdxList');
[~,ind] = max([s.Area]);
mask = false(size(ima));
mask(s(ind).PixelIdxList) = true;
