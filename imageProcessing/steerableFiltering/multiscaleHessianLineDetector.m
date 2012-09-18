%[response, theta, nms, scaleindex] = multiscaleHessianLineDetector(input, sigmaVect)

% Francois Aguet, Oct. 13, 2011

function [response, theta, nms, scaleindex] = multiscaleHessianLineDetector(input, sigmaVect)

[ny,nx] = size(input);
ns = numel(sigmaVect);

response = cell(1,ns);
theta = cell(1,ns);

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
    theta_s = zeros(ny,nx);
    response_s = zeros(ny,nx);
    for i = 1:nx*ny
        H = [f_xx(i) f_xy(i);
            f_xy(i) f_yy(i)];
        [V,D] = eig(-H); % peak of 2nd derivative is negative
        eigenValues = diag(D);
        maxIdx = find(eigenValues==max(eigenValues), 1, 'first');
        theta_s(i) = atan(V(2,maxIdx)/V(1,maxIdx));
        response_s(i) = eigenValues(maxIdx);
    end
    theta{si} = theta_s;
    response{si} = s^2 * response_s;
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
