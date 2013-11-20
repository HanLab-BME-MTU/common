%[response, theta, nms, scaleindex] = multiscaleHessianLineDetector(input, sigmaVect)

% Francois Aguet, Oct. 13, 2011

%Revised Hunter Elliott, Nov 2013;

function [response, theta, nms, scaleindex] = multiscaleHessianLineDetector(input, sigmaVect)

[ny,nx] = size(input);
ns = numel(sigmaVect);

response = cell(1,ns);
theta = cell(1,ns);

eigVal = cell(1,ns);
eigVec = cell(1,ns);

for si = 1:ns
    s = sigmaVect(si);

    w = ceil(4*s);
    x = -w:w;
    
    % 1-D components required for filtering
    g = exp(-x.^2/(2*s^2)) / (sqrt(2*pi)*s);
    gx = -x/s^2 .* g;
    gxx = x.^2 .* g / s^4; % -1/s^2 term subtracted below
    
    % compute 3 basis templates
    inputXT = padarray(input, [w w], 'symmetric');
    f_blur = conv2(g, g, inputXT, 'valid') / s^2; % col, row kernel
    f_xx = conv2(g, gxx, inputXT, 'valid') - f_blur;
    f_xy = conv2(gx, gx, inputXT, 'valid');
    f_yy = conv2(gxx, g, inputXT, 'valid') - f_blur;
    
    % eigenvalues -> response
    eigVal{si} = zeros(ny,nx,2);
    eigVec{si} = zeros(ny,nx,2,2);
    %Quadratic solution to eigenvalue problem for 2x2 symmetric matrix of H:
    % lambda1/2 = (f_xx + f_yy +/- sqrt((f_xx - f_yy) .^2 + 4*f_xy .^2)) ./ 2;
    %             ^----------^     ^------------------------------------^
    %                Alpha                   beta
    alpha = (f_xx + f_yy)/2;    
    beta = sqrt((f_xx - f_yy) .^2 + 4*f_xy .^2)/2; 
    eigVal{si}(:,:,1) = -alpha - beta;%Flip sign because we want eigenvalues of -H
    eigVal{si}(:,:,2) = -alpha + beta;
    
    %Get non-unit eigenvectors - we only use the direction    
    eigVec{si}(:,:,1,1) = eigVal{si}(:,:,1) + f_yy;
    eigVec{si}(:,:,2,1) = -f_xy;
            
    eigVec{si}(:,:,1,2) = eigVal{si}(:,:,2) + f_yy;
    eigVec{si}(:,:,2,2) = -f_xy;
            
    %Second eigenvalue/vector will always be largest.
    response{si} = eigVal{si}(:,:,2) * s^2;
    theta{si} = atan(eigVec{si}(:,:,2,2) ./ eigVec{si}(:,:,1,2));        
    
end

maxResponse = response{1};
maxTheta = theta{1};
scaleindex = ones(ny,nx);
for si = 2:ns
    idx = response{si} > maxResponse;
    maxResponse(idx) = response{si}(idx);
    maxTheta(idx) = theta{si}(idx);
    scaleindex(idx) = si;
end

response = maxResponse;
theta = maxTheta;
nms = nonMaximumSuppression(response, theta);
