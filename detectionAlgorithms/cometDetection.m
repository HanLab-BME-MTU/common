% featuresInfo = cometDetection(img, mask, sigma, mode)
%
% Inputs :      img : input image
%              mask : cell mask
%             sigma : standard deviation of the Gaussian PSF
%            {mode} : parameters to estimate, default 'xyArtc'
%           {alpha} : 
%
% Outputs:  featuresInfo : output structure with anisotropic Gaussian
%                          parameters, standard deviations (compatible with
%                          Khuloud's tracker.
%
% Sylvain Berlemont, April 2011

function featuresInfo = cometDetection(img, mask, sigma, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('mask', @islogical);
ip.addRequired('sigma', @isscalar);
ip.addParamValue('mode', 'xyArtc', @ischar);
ip.addParamValue('alpha', 0.05, @isscalar);
ip.addParamValue('kSigma', 3, @isscalar);
ip.addParamValue('minDist', .25, @isscalar);

ip.parse(img, mask, sigma, varargin{:});
mode = ip.Results.mode;
alpha = ip.Results.alpha;
kSigma = ip.Results.kSigma;
minDist = ip.Results.minDist;

img = double(img);
[nrows ncols] = size(img);

% Filter image with laplacian
bandPassIso = filterLoG(img,sigma);
bandPassIso(bandPassIso < 0) = 0;
bandPassIso(~mask) = 0;

% Filter image with steerable filter
[~,T] = steerableFiltering(img,2,sigma);
 
% Compute the local maxima of the bandpass filtered images
locMaxIso = locmax2d(bandPassIso, [5 5]);
 
bw = blobSegmentThreshold(bandPassIso,0,0,mask);

locMaxIso(~bw) = 0;

indMax = find(locMaxIso);
[y x] = ind2sub(size(img), indMax);

P = zeros(size(y, 1), 7);
P(:,1) = x;
P(:,2) = y;
P(:,3) = img(indMax);
P(:,4) = sigma;     % sigmaX
P(:,5) = sigma; % sigmaY
P(:,6) = T(indMax);

% Subresolution detection
hside = ceil(kSigma * sigma);
npx = (2 * hside + 1)^2;
xmin = x - hside;
xmax = x + hside;
ymin = y - hside;
ymax = y + hside;

isValid = find(xmin >= 1 & xmax <= ncols & ymin >= 1 & ymax <= nrows);

xmin = xmin(isValid);
xmax = xmax(isValid);
ymin = ymin(isValid);
ymax = ymax(isValid);
P = P(isValid,:);

stdP = zeros(size(P));
stdR = zeros(size(P,1),1);

kLevel = norminv(1 - alpha / 2.0, 0, 1); % ~2 std above background

success = false(numel(xmin),1);

for iFeature = 1:numel(xmin)
        
    crop = img(ymin(iFeature):ymax(iFeature), xmin(iFeature):xmax(iFeature));
    P(iFeature,7) = min(crop(:)); % background
    P(iFeature,3) = P(iFeature,3) - P(iFeature,7); % amplitude above background
        
    [params, stdParams, ~, res] = fitAnisoGaussian2D(crop, ...
        [0, 0, P(iFeature,3), 3 * P(iFeature,4), P(iFeature,5), ...
        P(iFeature,6), P(iFeature,7)], mode);
        
    % TEST: position must remain in a confined area
    isValid = max(abs(params(1:2))) < hside;
    
    % TEST: sigmaX > 1
    isValid = isValid & params(4) > 1;
    
    % TEST: goodness-of-fit
    stdRes = std(res(:));
    [~, pval] = kstest(res(:) ./ stdRes, [], alpha);
    isValid = isValid & pval > alpha;

    % TEST: amplitude
    SE_sigma_r = (stdRes / sqrt(2*(npx-1))) * kLevel;
    sigma_A = stdParams(3);
    A_est = params(3);
    df2 = (npx - 1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
    scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2) / npx);
    T = (A_est - stdRes * kLevel) ./ scomb;    
    isValid = isValid & (1 - tcdf(T, df2)) < alpha;
  
    % TEST: extreme value of 
    isValid = isValid & params(4) < 10 * sigma;
    
    success(iFeature) = isValid;
    
    P(iFeature,1) = P(iFeature,1) + params(1);
    P(iFeature,2) = P(iFeature,2) + params(2);
    P(iFeature,3) = params(3);
    P(iFeature,4) = params(4);
    P(iFeature,5) = params(5);
    P(iFeature,6) = params(6);
    P(iFeature,7) = params(7);
    
    stdP(iFeature,1) = stdParams(1);
    stdP(iFeature,2) = stdParams(2);
    stdP(iFeature,3) = stdParams(3);
    stdP(iFeature,4) = stdParams(4);
    stdP(iFeature,6) = stdParams(5);
    stdP(iFeature,7) = stdParams(6);
    
    stdR(iFeature) = stdRes;
end

P = P(success,:);
stdP = stdP(success,:);

% Remove any detection which has been localised at the same position
isValid = true(size(P,1),1);
idxKD = KDTreeBallQuery(P(:,1:2), P(:,1:2), repmat(minDist, size(P,1), 1));
idxKD = idxKD(cellfun(@(x) length(x)>1, idxKD));
    
for k = 1:length(idxKD);
    stdRes = stdR(idxKD{k});
    isValid(idxKD{k}(stdRes ~= min(stdRes))) = false;
end

P = P(isValid,:);
stdP = stdP(isValid,:);

featuresInfo.xCoord = [P(:,1), stdP(:,1)];
featuresInfo.yCoord = [P(:,2), stdP(:,2)];
featuresInfo.amp = [P(:,3), stdP(:,3)];
featuresInfo.sigmaX = [P(:,4), stdP(:,4)];
featuresInfo.sigmaY = [P(:,5), stdP(:,5)];
featuresInfo.theta = [P(:,6), stdP(:,6)];
featuresInfo.bkg = [P(:,7), stdP(:,7)];
