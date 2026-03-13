% pStruct = fitGaussianMixtures2D(img, x, y, A, sigma, c, mode, varargin)
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
% Output: pStruct: structure with fields:
%                  x : estimated x-positions
%                  y : estimated y-positions
%                  A : estimated amplitudes
%                  s : estimated standard deviations of the PSF
%                  c : estimated background intensities
%             x_pstd : standard deviations, estimated by error propagation
%             y_pstd :
%             A_pstd :
%             s_pstd :
%             c_pstd :
%            sigma_r : standard deviation of the background (residual)
%         SE_sigma_r : standard error of sigma_r
%            pval_Ar : p-value of amplitude vs. background noise test
%
% Francois Aguet, March 28 2011 (last modified: Feb 5 2013)
% parfor conversion: Sangyoon Han, 2026

function pStruct = fitGaussianMixtures2D(img, x, y, A, sigma, c, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('x');
ip.addRequired('y');
ip.addRequired('A');
ip.addRequired('sigma');
ip.addRequired('c');
ip.addParamValue('Alpha',  0.05, @isscalar);
ip.addParamValue('AlphaT', 0.05, @isscalar);
ip.addParamValue('Mask',   [],   @islogical);
ip.addParamValue('maxM',   5,    @isscalar);
ip.parse(img, x, y, A, sigma, c, varargin{:});

np    = length(x);
alpha = ip.Results.Alpha;
alphaT= ip.Results.AlphaT;
maxM  = ip.Results.maxM;

sigma = ip.Results.sigma;
if numel(sigma) == 1
    sigma = sigma * ones(1, np);
end

if ~isempty(ip.Results.Mask)
    labels = bwlabel(ip.Results.Mask);
else
    labels = zeros(size(img));
end

% --- Pre-allocate output cell arrays (parfor-safe) -------------------------
x_out      = cell(1, np);   y_out      = cell(1, np);
A_out      = cell(1, np);   s_out      = cell(1, np);   c_out      = cell(1, np);
xp_out     = cell(1, np);   yp_out     = cell(1, np);
Ap_out     = cell(1, np);   sp_out     = cell(1, np);   cp_out     = cell(1, np);
xi_out     = cell(1, np);   yi_out     = cell(1, np);
sigr_out   = cell(1, np);   SEsigr_out = cell(1, np);
RSS_out    = cell(1, np);
pval_out   = cell(1, np);   hval_AD_out= cell(1, np);   hval_Ar_out= cell(1, np);
% mixtureIndex: store p (1-based) for mixtures, 0 for singles, -1 for empty
mixIdx_out = cell(1, np);

[ny, nx_img] = size(img);
kLevel = norminv(1 - alpha/2.0, 0, 1);
iRange = [min(img(:)) max(img(:))];

sigma_max = max(sigma);
w2 = ceil(2*sigma_max);
w3 = ceil(3*sigma_max);
w4 = ceil(4*sigma_max);

% annular mask for background estimation
[xm, ym] = meshgrid(-w4:w4);
r = sqrt(xm.^2 + ym.^2);
annularMask = double(r <= w4 & r >= w3);

xi = round(x);
yi = round(y);

% Snapshot loop variables for parfor broadcast
A_in = A;
c_in = c;
x_in = x;
y_in = y;
sigma_in = sigma;

parfor p = 1:np %#ok<*PFBNS>

    % Skip border points
    if xi(p) <= w4 || xi(p) > nx_img-w4 || yi(p) <= w4 || yi(p) > ny-w4
        continue
    end

    % label window ? mask out other nearby spots
    maskWindow = labels(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
    maskWindow(maskWindow == maskWindow(w4+1, w4+1)) = 0;

    % background estimation
    cmask = annularMask;
    cmask(maskWindow ~= 0) = 0;
    window = img(yi(p)-w4:yi(p)+w4, xi(p)-w4:xi(p)+w4);
    if isempty(c_in)
        c0 = mean(window(cmask == 1));
    else
        c0 = c_in(p);
    end

    % mask out other components
    window(maskWindow ~= 0) = NaN;
    npx = sum(isfinite(window(:)));

    if isempty(A_in)
        A0 = max(window(:)) - c0;
    else
        A0 = A_in(p);
    end

    % --- Initial single-Gaussian fit (reduced model) -----------------------
    [prm_f, prmStd_f, ~, res_f] = fitGaussian2D(window, ...
        [x_in(p)-xi(p)  y_in(p)-yi(p)  A0  sigma_in(p)  c0], 'xyAc');
    prmStd_f = [prmStd_f(1:3) 0 prmStd_f(4)]; % insert zero for fixed sigma std

    RSS_r  = res_f.RSS;
    p_r    = 4;   % #params in reduced model (x,y,A,c)
    i      = 1;   % mixture iteration counter

    pval       = 0;
    validBounds= true;

    % --- Iteratively add Gaussian components --------------------------------
    while i < maxM && pval < 0.05 && validBounds

        i = i + 1;
        prm_r    = prm_f;
        prmStd_r = prmStd_f;
        res_r    = res_f;

        % New component: seed at max residual
        [y0, x0] = find(res_r.data == max(res_r.data(:)), 1, 'first');
        [prm_f, prmStd_f, ~, res_f] = fitGaussianMixture2D(window, ...
            [x0-w4-1  y0-w4-1  max(res_r.data(:))  prm_r], 'xyAc');

        RSS_f = res_f.RSS;
        p_f   = p_r + 3;   % 3 new params (x,y,A) per component

        % F-test
        T    = (RSS_r - RSS_f) / RSS_f * (npx - p_f - 1) / (p_f - p_r);
        pval = 1 - fcdf(T, p_f-p_r, npx-p_f-1);

        p_r   = p_f;
        RSS_r = RSS_f;

        % Restrict radius ? ignore out-of-window components
        if min(prm_f(1:3:end-2)) < -w2 || max(prm_f(1:3:end-2)) > w2 || ...
           min(prm_f(2:3:end-2)) < -w2 || max(prm_f(2:3:end-2)) > w2
            validBounds = false;
        end
    end

    ng    = i - 1;  % # Gaussians in final model
    x_est = prm_r(1:3:end-2);
    y_est = prm_r(2:3:end-2);
    A_est = prm_r(3:3:end-2);

    % Accept fit if multi-component OR single within bounds
    if ng > 1 || (x_est > -w2 && x_est < w2 && ...
                  y_est > -w2 && y_est < w2 && A_est < 2*diff(iRange))

        x_out{p}   = xi(p) + x_est;
        y_out{p}   = yi(p) + y_est;
        A_out{p}   = A_est;
        s_out{p}   = repmat(prm_r(end-1), [1 ng]);
        c_out{p}   = repmat(prm_r(end),   [1 ng]);

        xp_out{p}  = prmStd_r(1:3:end-2);
        yp_out{p}  = prmStd_r(2:3:end-2);
        Ap_out{p}  = prmStd_r(3:3:end-2);
        sp_out{p}  = repmat(prmStd_r(end-1), [1 ng]);
        cp_out{p}  = repmat(prmStd_r(end),   [1 ng]);

        xi_out{p}  = repmat(xi(p), [1 ng]);
        yi_out{p}  = repmat(yi(p), [1 ng]);

        sigr_out{p}  = repmat(res_r.std, [1 ng]);
        RSS_out{p}   = repmat(res_r.RSS, [1 ng]);

        SE_sigma_r         = res_r.std / sqrt(2*(npx-1));
        SEsigr_out{p}      = repmat(SE_sigma_r, [1 ng]);
        SE_sigma_r_kLevel  = SE_sigma_r * kLevel;

        hval_AD_out{p} = repmat(res_r.hAD, [1 ng]);

        % Amplitude significance test (one per mixture component)
        pval_p   = zeros(1, ng);
        hval_p   = false(1, ng);
        for ig = 1:ng
            sigma_A = Ap_out{p}(ig);
            if ~isfinite(sigma_A) || sigma_A < eps
                sigma_A = eps;
            end
            Aig  = A_out{p}(ig);
            df2  = (npx-1) * (sigma_A^2 + SE_sigma_r_kLevel^2)^2 / ...
                             (sigma_A^4  + SE_sigma_r_kLevel^4);
            scomb = sqrt((sigma_A^2 + SE_sigma_r_kLevel^2) / npx);
            T_ig  = (Aig - res_r.std*kLevel) / scomb;
            pval_p(ig)  = tcdf(-T_ig, df2);
            hval_p(ig)  = pval_p(ig) < alphaT;
        end
        pval_out{p}    = pval_p;
        hval_Ar_out{p} = hval_p;

        % Mixture index: use p as a unique placeholder for multi-component
        % spots; re-number sequentially after the loop.
        if ng > 1
            mixIdx_out{p} = p * ones(1, ng);   % unique per worker slot
        else
            mixIdx_out{p} = zeros(1, ng);       % 0 = single Gaussian
        end
    end
end % parfor

% --- Re-number mixtureIndex sequentially (1, 2, 3, ...) -------------------
% parfor used p as a placeholder; now compress to consecutive integers.
counter = 0;
for p = 1:np
    if ~isempty(mixIdx_out{p}) && any(mixIdx_out{p} > 0)
        counter = counter + 1;
        mixIdx_out{p}(:) = counter;
    end
end

% --- Assemble pStruct -------------------------------------------------------
pStruct = struct( ...
    'x',          [x_out{:}], ...
    'y',          [y_out{:}], ...
    'A',          [A_out{:}], ...
    's',          [s_out{:}], ...
    'c',          [c_out{:}], ...
    'x_pstd',     [xp_out{:}], ...
    'y_pstd',     [yp_out{:}], ...
    'A_pstd',     [Ap_out{:}], ...
    's_pstd',     [sp_out{:}], ...
    'c_pstd',     [cp_out{:}], ...
    'x_init',     [xi_out{:}], ...
    'y_init',     [yi_out{:}], ...
    'sigma_r',    [sigr_out{:}], ...
    'SE_sigma_r', [SEsigr_out{:}], ...
    'RSS',        [RSS_out{:}], ...
    'pval_Ar',    [pval_out{:}], ...
    'hval_AD',    [hval_AD_out{:}], ...
    'hval_Ar',    [hval_Ar_out{:}], ...
    'mixtureIndex',[mixIdx_out{:}]);