function ima = ib3spline2D(coeffs, imaSize)

if nargin > 1 && ~isempty(imaSize)
    [ny,nx] = size(coeffs);
else
    ny = imaSize(1);
    nx = imaSize(2);
end

[nyCoeffs, nxCoeffs] = size(coeffs);

[X,Y] = meshgrid(1:nx,1:ny);

% Rescale coordinates
xScaled = X / (nx / nxCoeffs);
yScaled = Y / (ny / nyCoeffs);

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
coeffsPadded = padarrayXT(coeffs, [2 2], 'symmetric', 'both');

% Compute interpolation for each X,Y = 
% sum_k=1^4 yWeight(k) sum_l=1^4 xWeight(l) BSpadded(yIndex(l), xIndex(l))

ima = zeros(ny,nx);

tmp = zeros(ny,nx);

for k = 1:4
    tmp(:) = 0;
    
    yShifted = yIndex(:,:,k) + 2;
    
    for l = 1:4
        xShifted = xIndex(:,:,l) + 2;
        
        ind = sub2ind([nyCoeffs, nxCoeffs], yShifted(:), xShifted(:));
        
        tmp = tmp + xWeight(:,:,l) .* reshape(coeffsPadded(ind), [ny, nx]);
        
    end
    
    ima = ima + yWeight(:,:,k) .* tmp;
end