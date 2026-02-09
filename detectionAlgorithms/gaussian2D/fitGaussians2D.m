% pStruct = fitGaussians2D(img, x, y, A, sigma, c, mode, varargin)
%
% Inputs:   img : input image
%             x : initial (or fixed) x-positions
%             y : initial (or fixed) y-positions
%             A : initial (or fixed) amplitudes
%             s : initial (or fixed) Gaussian PSF standard deviations
%             c : initial (or fixed) background intensities
%
% Options:
%          mode : selector for optimization parameters, any of 'xyAsc'. Default: 'xyAc'
%
% Options ('specifier', value):
%          'Mask' : mask of spot locations
%    'ConfRadius' : Confinement radius for valid fits. Default: ceil(2*sigma)
%    'WindowSize' : Size of the support used for the fit, specified as half-width;
%                   i.e., for a window of 15x15, enter 7. Default: ceil(4*sigma)
%
% Output: pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
%                  c : estimated background intensities
%
%             x_pstd : standard deviations, estimated by error propagation
%             y_pstd : "
%             A_pstd : "
%             s_pstd : "
%             c_pstd : "
%            sigma_r : standard deviation of the background (residual)
%         SE_sigma_r : standard error of sigma_r
%            pval_Ar : p-value of an amplitude vs. background noise test (p < 0.05 -> significant amplitude)
%
%
% Usage for a single-channel image with mask and fixed sigma:
% fitGaussians2D(img, x, y, A, sigma, c, 'xyAc', 'mask', mask);

% Francois Aguet, March 28 2011 (last updated: Sep 30 2013)

function pStruct = fitGaussians2D(img, x, y, A, sigma, c, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('x');
ip.addRequired('y');
ip.addRequired('A');
ip.addRequired('sigma');
ip.addRequired('c');
ip.addOptional('mode', 'xyAc', @ischar);
ip.addParamValue('Alpha', 0.05, @isscalar);
ip.addParamValue('AlphaT', 0.05, @isscalar);
ip.addParamValue('Mask', [], @islogical);
ip.addParamValue('ConfRadius', []);
ip.addParamValue('WindowSize', []);
ip.parse(img, x, y, A, sigma, c, varargin{:});

np = length(x);
sigma = ip.Results.sigma;
if numel(sigma)==1
    sigma = sigma*ones(1,np);
end
mode = ip.Results.mode;
if ~isempty(ip.Results.Mask)
    labels = bwlabel(ip.Results.Mask);
else
    labels = zeros(size(img));
end

pStruct = struct('x', [], 'y', [], 'A', [], 's', [], 'c', [],...
    'x_pstd', [], 'y_pstd', [], 'A_pstd', [], 's_pstd', [], 'c_pstd', [],...
    'x_init', [], 'y_init', [],...
    'sigma_r', [], 'SE_sigma_r', [], 'RSS', [], 'pval_Ar', [], 'mask_Ar', [], 'hval_Ar', [], 'hval_AD', []);

xi = round(x);
yi = round(y);
[ny,nx] = size(img);

kLevel = norminv(1-ip.Results.Alpha/2.0, 0, 1); % ~2 std above background

iRange = [min(img(:)) max(img(:))];

estIdx = regexpi('xyAsc', ['[' mode ']']);


% initialize pStruct arrays
pStruct.x = NaN(1,np);
pStruct.y = NaN(1,np);
pStruct.A = NaN(1,np);
pStruct.s = NaN(1,np);
pStruct.c = NaN(1,np);

pStruct.x_pstd = NaN(1,np);
pStruct.y_pstd = NaN(1,np);
pStruct.A_pstd = NaN(1,np);
pStruct.s_pstd = NaN(1,np);
pStruct.c_pstd = NaN(1,np);

pStruct.x_init = reshape(xi, [1 np]);
pStruct.y_init = reshape(yi, [1 np]);

pStruct.sigma_r = NaN(1,np);
pStruct.SE_sigma_r = NaN(1,np);
pStruct.RSS = NaN(1,np);

pStruct.pval_Ar = NaN(1,np);
pStruct.mask_Ar = zeros(1,np);

pStruct.hval_AD = false(1,np);
pStruct.hval_Ar = false(1,np);


sigma_max = max(sigma);
w2 = ip.Results.ConfRadius;
if isempty(w2)
    w2 = ceil(2*sigma_max);
end
w4 = ip.Results.WindowSize;
if isempty(w4)
    w4 = ceil(4*sigma_max);
end

% for background estimation
if isempty(c)
    % mask template: ring with inner radius w3, outer radius w4
    [xm,ym] = meshgrid(-w4:w4);
    r = sqrt(xm.^2+ym.^2);
    annularMask = zeros(size(r));
    annularMask(r<=ceil(4*sigma_max) & r>=ceil(3*sigma_max)) = 1;
end

g = exp(-(-w4:w4).^2/(2*sigma_max^2));
g = g'*g;
g = g(:);

T = zeros(1,np);
df2 = zeros(1,np);
for p = 1:np
    
    % ignore points in border
    if (xi(p)>w4 && xi(p)<=nx-w4 && yi(p)>w4 && yi(p)<=ny-w4)
        
        % label mask
        maskWindow = labels(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
        maskWindow(maskWindow==maskWindow(w4+1,w4+1)) = 0;
        
        window = img(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);

        % estimate background        
        if isempty(c)
            cmask = annularMask;
            cmask(maskWindow~=0) = 0;
            c_init = mean(window(cmask==1));
        else
            c_init = c(p);
        end
        
        % set any other components to NaN
        window(maskWindow~=0) = NaN;
        npx = sum(isfinite(window(:)));
        
        if npx >= 10 % only perform fit if window contains sufficient data points
        
            if isa(window,'integer')
                window  = double(window);
            end
            
            % fit
            if isempty(A)
                A_init = max(window(:))-c_init;
            else
                A_init = A(p);
            end
            
            [prm, prmStd, ~, res, ~] = fitGaussian2D(window, [x(p)-xi(p) y(p)-yi(p) A_init sigma(p) c_init], mode);
            % ---- Begin fallback for A-std if covariance is singular ----
            % Only attempt if A is estimated by 'mode'
            if any(mode == 'A')
                if numel(prmStd) >= 3 && (isnan(prmStd(3)) || ~isfinite(prmStd(3)) || prmStd(3) <= 0)
                    % Build the Gaussian shape on the same window the fit used
                    [H, W] = size(window);
                    cx = (W+1)/2; cy = (H+1)/2;                 % window center (origin used by the wrapper)
                    [Xg, Yg] = meshgrid((1:W)-cx, (1:H)-cy);
                    x0 = prm(1); y0 = prm(2); s = prm(4);
            
                    % Valid pixels (same criterion as fit used)
                    mask = isfinite(window);
                    if ~any(mask(:))
                        % nothing we can do; leave NaN
                    else
                        phi = exp(-((Xg - x0).^2 + (Yg - y0).^2) / (2*s^2));
                        phi = phi(mask);
            
                        % m = number of residuals actually used; p_est = #parameters estimated
                        m = sum(mask(:));
                        p_est = numel(intersect(lower(mode), 'xyasc')); % estimated params count
                        dof = max(m - p_est, 1);                         % protect against 0/neg
            
                        % sigma_r^2 from residuals already returned by the MEX
                        % res.RSS should also be available; prefer it if present
                        if isfield(res, 'RSS') && isfinite(res.RSS)
                            sigma2 = res.RSS / dof;
                        else
                            sigma2 = res.std.^2; % close enough if RSS not exposed
                        end
            
                        denom = sum(phi.^2);
                        if denom > eps
                            stdA = sqrt(sigma2 / denom);
                            prmStd(3) = stdA;  % overwrite NaN with fallback
                        end
                    end
                end
            end
            % ---- End fallback for A-std ----
            
            dx = prm(1);
            dy = prm(2);

            % --- BEGIN: compatibility tests expected by upstream code ---
            
            r = res.data(:);
            
            % 1) Anderson–Darling normality test (preferred by your pipeline)
            if exist('adtest','file') == 2
                try
                    warning('off','stats:adtest:OutOfRangePLow');
                    warning('off','stats:adtest:OutOfRangePHigh');

                    [hAD, pAD] = adtest(r);

                    warning('on','stats:adtest:OutOfRangePLow');
                    warning('on','stats:adtest:OutOfRangePHigh');
                catch
                    hAD = NaN; pAD = NaN;
                end
            else
                % Fallback if adtest is unavailable
                hAD = NaN; pAD = NaN;
            end
            res.hAD = hAD;
            res.pAD = pAD;
            
            % 2) (Optional) Also provide KS flags some codebases expect
            if exist('kstest','file') == 2
                try
                    % Test against N(mean, std^2) inferred from residuals
                    rZ = (r - mean(r,'omitnan')) / std(r,[],'omitnan');
                    [hKS, pKS] = kstest(rZ);
                catch
                    hKS = NaN; pKS = NaN;
                end
            else
                hKS = NaN; pKS = NaN;
            end
            res.hKS = hKS;
            res.pKS = pKS;
            
            % --- END: compatibility tests ---

            
            % exclude points where localization failed
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
                
                % --- BEGIN robust T/df2 computation ---
                
                % Guard residual std (recompute from residuals if NaN or empty)
                if ~isfinite(pStruct.sigma_r(p))
                    rtmp = res.data(:);
                    rtmp = rtmp(isfinite(rtmp));
                    if ~isempty(rtmp)
                        pStruct.sigma_r(p) = std(rtmp, 0); % unbiased by default (n-1)
                    else
                        pStruct.sigma_r(p) = NaN;
                    end
                end
                
                % Standard error of sigma_r (guard npx and finiteness)
                if isfinite(pStruct.sigma_r(p)) && npx > 1
                    pStruct.SE_sigma_r(p) = pStruct.sigma_r(p) / sqrt(max(2*(npx-1),1));
                else
                    pStruct.SE_sigma_r(p) = NaN;
                end
                
                SE_sigma_r = pStruct.SE_sigma_r(p) * kLevel;
                
                % Pull sigma_A (std of amplitude estimate) from stdVect; if A not estimated, it is 0
                sigma_A = stdVect(3);
                
                % If A was fixed or prmStd missing, avoid exact zero to keep formulas defined
                if ~isfinite(sigma_A) || sigma_A < eps
                    % If A is fixed, the only uncertainty is from background; treat sigma_A ~ 0
                    % but keep it tiny to avoid 0/0 in Welch df:
                    sigma_A = eps;
                end
                
                A_est = prm(3);
                
                % Welch–Satterthwaite df for sum of two variances:
                num = (sigma_A.^2 + SE_sigma_r.^2).^2;
                den = (sigma_A.^4 + SE_sigma_r.^4);
                
                if isfinite(num) && isfinite(den) && den > 0 && npx > 1
                    df2(p) = (npx - 1) * num / den;
                else
                    df2(p) = NaN;
                end
                
                % Combined standard error for 1-sample t on (A_est vs k*sigma_r)
                scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2) / max(npx,1));
                
                if isfinite(scomb) && scomb > 0 && isfinite(pStruct.sigma_r(p))
                    T(p) = (A_est - pStruct.sigma_r(p)*kLevel) / scomb;
                else
                    T(p) = NaN;
                end
                
                % Mask count used upstream (keep your original logic; guard res.std)
                if isfinite(pStruct.sigma_r(p))
                    pStruct.mask_Ar(p) = sum(A_est * g > pStruct.sigma_r(p) * kLevel);
                else
                    pStruct.mask_Ar(p) = 0;
                end
                
                % Keep the AD flag you already stored
                pStruct.hval_AD(p) = res.hAD;
                
                % --- END robust T/df2 computation ---                
                pStruct.SE_sigma_r(p) = res.std/sqrt(2*(npx-1));
                SE_sigma_r = pStruct.SE_sigma_r(p) * kLevel;
                
                pStruct.hval_AD(p) = res.hAD;
                
                % H0: A <= k*sigma_r
                % H1: A > k*sigma_r
                sigma_A = stdVect(3);
                A_est = prm(3);
                df2(p) = (npx-1) * (sigma_A.^2 + SE_sigma_r.^2).^2 ./ (sigma_A.^4 + SE_sigma_r.^4);
                scomb = sqrt((sigma_A.^2 + SE_sigma_r.^2)/npx);
                T(p) = (A_est - res.std*kLevel) ./ scomb;
                pStruct.mask_Ar(p) = sum(A_est*g>res.std*kLevel);
            end
        end
    end
end
% 1-sided t-test: A_est must be greater than k*sigma_r
pStruct.pval_Ar = tcdf(-T, df2);
pStruct.hval_Ar = pStruct.pval_Ar < ip.Results.AlphaT;
