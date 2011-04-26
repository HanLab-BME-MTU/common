% pStruct = fitGaussians2D(img, x, y, A, sigma, c, mode, varargin)
%
% Inputs:   img : input image
%             x : initial (or fixed) x-positions
%             y : initial (or fixed) y-positions
%             A : initial (or fixed) amplitudes
%             s : initial (or fixed) Gaussian PSF standard deviations
%             c : initial (or fixed) background intensities
%          mode : string selector for optimization parameters, any of 'xyAsc'
%
% Optional inputs : ('Mask', mask) pair with a mask of spot locations
%
%
% Output: pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
%                  c : estimated background intensities
%
%             x_pstd : standard deviations, estimated by error propagation
%             y_pstd :
%             A_pstd :
%             s_pstd :
%             c_pstd :
%            sigma_r : standard deviation of the background (residual)
%         SE_sigma_r : standard error of sigma_r
%            pval_KS : p-value of the Kolomogorv-Smirnov test on the residuals (p > 0.05 -> Gaussian)
%            pval_Ar : p-value of an amplitude vs. background noise test (p > 0.95 -> significant amplitude)
%
%
% Usage for a single-channel img with mask and fixed sigma:
% fitGaussians2D(img, x_v, y_v, 'sigma', sigma_v, 'mask', mask);

% Francois Aguet, March 28 2011 (last modified: April 5 2011)

function pStruct = fitGaussians2D(img, x, y, A, sigma, c, mode, varargin)

np = length(x);

pStruct = struct('x', [], 'y', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 's_pstd', [], 'c_pstd', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_KS', [], 'pval_Ar', []);

xi = round(x);
yi = round(y);
[ny,nx] = size(img);

if mod(length(varargin),2)~=0
    error('Optional arguments need to be entered as pairs.');
end

idx = find(strcmpi(varargin, 'alpha'));
if ~isempty(idx)
    alpha = varargin{idx+1};
else
    alpha = 0.05;
end

idx = find(strcmpi(varargin, 'mask'));
if ~isempty(idx)
    mask = varargin{idx+1};
    labels = bwlabel(mask);
else
    labels = zeros(size(img));
end

kLevel = norminv(1-alpha/2.0, 0, 1); % ~2 std above background

iRange = [squeeze(min(min(img))) squeeze(max(max(img)))];

estIdx = regexpi('xyAsc', ['[' mode ']']);


% initialize pStruct arrays
pStruct.x = NaN(1,np);
pStruct.y = NaN(1,np);
pStruct.A = NaN(1,np);


pStruct.x_pstd = NaN(1,np);
pStruct.y_pstd = NaN(1,np);
pStruct.A_pstd = NaN(1,np);
% if fitSigma
pStruct.s = NaN(1,np);
pStruct.s_pstd = NaN(1,np);
% end
pStruct.c = NaN(1,np);
pStruct.c_pstd = NaN(1,np);
pStruct.sigma_r = NaN(1,np);
pStruct.SE_sigma_r = NaN(1,np);
pStruct.RSS = NaN(1,np);

pStruct.pval_Ar = NaN(1,np);
pStruct.pval_KS = NaN(1,np);

sigma_max = max(sigma);
w2 = ceil(2*sigma_max);
w3 = ceil(3*sigma_max);
w4 = ceil(4*sigma_max);

% mask template: ring with inner radius w3, outer radius w4
[xm,ym] = meshgrid(-w4:w4);
r = sqrt(xm.^2+ym.^2);
annularMask = zeros(size(r));
annularMask(r<=w4 & r>=w3) = 1;


% indexes for cross-correlation coefficients
% n = length(mode);
% idx = pcombs(1:n);
% i = idx(:,1);
% j = idx(:,2);
% ij = i+n*(j-1);
% ii = i+n*(i-1);
% jj = j+n*(j-1);


for p = 1:np
    
    % ignore points in border
    if (xi(p)>w4 && xi(p)<=nx-w4 && yi(p)>w4 && yi(p)<=ny-w4)
        
        % label mask
        maskWindow = labels(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
        
        % estimate background
        cmask = annularMask;
        cmask(maskWindow~=0) = 0;
        window = img(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        if isempty(c)
            c_init = mean(window(cmask==1));
        else
            c_init = c(p);
        end
        
        % set any other components to NaN
        window(maskWindow~=0) = NaN;
        
        % standard deviation of the background within annular mask
        %bgStd = nanstd(window(annularMask==1));
        
        npx = sum(isfinite(window(:)));
        
        % fit
        if isempty(A)
            A_init = max(window(:))-c_init;
        else
            A_init = A(p);
        end

        [prm, prmStd, ~, res] = fitGaussian2D(window, [x(p)-xi(p) y(p)-yi(p) A_init sigma(p) c_init], mode);
        %K = corrFromC(C,ij,ii,jj);
        
        dx = prm(1);
        dy = prm(2);
        
        % eliminate points where localization failed or which are close to image border
        if (dx > -w2 && dx < w2 && dy > -w2 && dy < w2 && prm(3)<2*diff(iRange))
            
            pStruct.x(p) = xi(p) + dx;
            pStruct.y(p) = yi(p) + dy;
            pStruct.A(p) = prm(3);
            pStruct.s(p) = prm(4);
            pStruct.c(p) = prm(5);
            
            stdVect = zeros(1,5);
            stdVect(estIdx) = prmStd;
            
            pStruct.x_pstd(p) = stdVect(1);
            pStruct.y_pstd(p) = stdVect(2);
            pStruct.A_pstd(p) = stdVect(3);
            pStruct.s_pstd(p) = stdVect(4);
            pStruct.c_pstd(p) = stdVect(5);
            
            pStruct.sigma_r(p) = res.std;
            pStruct.RSS(p) = res.RSS;
            
            pStruct.SE_sigma_r(p) = res.std/sqrt(2*(npx-1));
            SE_sigma_r = pStruct.SE_sigma_r(p) * kLevel;
            
            pStruct.pval_KS(p) = res.pval;
            
            sigma_A = stdVect(3);
            A_est = prm(3);
            df2 = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
            scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
            T = (A_est - res.std*kLevel) ./ scomb;
            pStruct.pval_Ar(p) = tcdf(T, df2);
        end
    end
end


% function K = corrFromC(C,ij,ii,jj)
% n = size(C,1);
% K = zeros(n,n);
%
% K(ij) = C(ij) ./ sqrt(C(ii).*C(jj));
% remaining components are redundant
% K = K + K';
% K(sub2ind([n n], 1:n, 1:n)) = 1;
