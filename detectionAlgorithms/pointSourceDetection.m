%[pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, sigma, mode)
%
% Inputs :      img : input image
%             sigma : standard deviation of the Gaussian PSF
%            {mode} : parameters to estimate, default 'xyAc'
%
% Outputs:  pstruct : output structure with Gaussian parameters, standard deviations, p-values
%              mask : mask of significant (in amplitude) pixels
%             imgLM : image of local maxima
%            imgLoG : Laplacian of Gaussian filtered image

% Francois Aguet, April 2-5 2011

function [pstruct, mask, imgLM, imgLoG] = pointSourceDetection(img, sigma, mode)

if nargin<3
    mode = 'xyAc';
end

% Gaussian kernel
w = ceil(4*sigma);
x = -w:w;
g = exp(-x.^2/(2*sigma^2));
u = ones(1,length(x));

imgXT = padarrayXT(img, [w w], 'symmetric');

% convolutions
fg = conv2(g', g, imgXT, 'valid');
fu = conv2(u', u, imgXT, 'valid');
fu2 = conv2(u', u, imgXT.^2, 'valid');

% Laplacian of Gaussian
gx2 = g.*x.^2;
imgLoG = 2*fg/sigma^2 - (conv2(g, gx2, imgXT, 'valid')+conv2(gx2, g, imgXT, 'valid'))/sigma^4;
imgLoG = imgLoG / (2*pi*sigma^2);

g = g'*g;
n = numel(g);

gsum = sum(g(:));
g2sum = sum(g(:).^2);

A_est = (fg - gsum*fu/n) / (g2sum - gsum^2/n);
c_est = (fu - A_est*gsum)/n;

J = [g(:) ones(n,1)]; % g_dA g_dc
C = inv(J'*J);

f_c = fu2 - 2*c_est.*fu + n*c_est.^2; % f-c
RSS = A_est.^2*g2sum - 2*A_est.*(fg - c_est*gsum) + f_c;
sigma_e2 = RSS/(n-3);

sigma_A = sqrt(sigma_e2*C(1,1));

% std(residuals)
sigma_res = sqrt((RSS - (A_est*gsum+n*c_est - fu)/n)/(n-1));

kLevel = norminv(0.95,0,1);

SE_sigma_c = sigma_res/sqrt(2*(n-1)) * kLevel;
df2 = (n-1) * (sigma_A.^2 + SE_sigma_c.^2).^2 ./ (sigma_A.^4 + SE_sigma_c.^4);
scomb = sqrt((sigma_A.^2 + SE_sigma_c.^2)/n);
T = (A_est - sigma_res*kLevel) ./ scomb;
pval = tcdf(real(T), df2);

% mask of admissible positions for local maxima
mask = pval > 0.95;
% mask = bwmorph(mask,'dilate');

% local maxima
imgLM = locmax2d(imgLoG, 5) .* mask;
[lmy, lmx] = find(imgLM~=0);
lmIdx = sub2ind(size(img), lmy, lmx);

% run localization on local maxima
pstruct = fitGaussians2D(img, lmx, lmy, A_est(lmIdx), sigma*ones(1,length(lmIdx)), c_est(lmIdx), mode);

% eliminate isignificant amplitudes
idx = [pstruct.pval_Ar] > 0.95;

fnames = fieldnames(pstruct);
for k = 1:length(fnames)
    pstruct.(fnames{k}) = pstruct.(fnames{k})(idx);
end
pstruct.isPSF = pstruct.pval_KS > 0.05;
