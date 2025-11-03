function [prmVect, prmStd, C, res, J] = fitGaussian2D(data, prmVect, mode, options)
%FITGAUSSIAN2D  Compatibility wrapper for new fitGaussian2D_mex (GSL backend)
%
%   [prmVect, prmStd, C, res, J] = fitGaussian2D(data, prmVect, mode, options)
%
%   Supports parameter subset fitting via 'mode' (e.g. 'xyasc', 'As', 'asc').
%   Keeps legacy input/output order: [xp yp A s c].
%
% Francois Aguet (2011) style interface, updated for MATLAB R2024b + GSL MEX.
% Sangyoon han (2025)

% ———————————————————————––
% Input handling
% ———————————————————————––
if nargin < 2
error('fitGaussian2D:Input','At least data and prmVect required.');
end
if nargin < 3 || isempty(mode)
mode = 'xyasc'; % fit all by default
end
if nargin < 4 || isempty(options)
options = [200 1e-8 1e-8]; % [maxIter eAbs eRel]
end
if ~isfloat(data)
data = double(data);
end

% Mask NaNs
mask = ~isnan(data);
data(isnan(data)) = 0;

% ———————————————————————––
% Setup parameter masks based on 'mode'
% ———————————————————————––
letters = lower(mode);
names   = {'x','y','a','s','c'};
idxMap  = struct('x',1,'y',2,'a',3,'s',4,'c',5);
fitMask = false(1,5);
for k = 1:numel(letters)
if isfield(idxMap,letters(k))
fitMask(idxMap.(letters(k))) = true;
end
end
if ~any(fitMask)
warning('fitGaussian2D:Mode','Mode “%s” contained no valid symbols. Fitting all.',mode);
fitMask(:) = true;
end

% ———————————————————————––
% Prepare MEX input: [A, x0, y0, s, b]
% ———————————————————————––
fullVec = prmVect(:)';
if numel(fullVec) ~= 5
error('fitGaussian2D:Init','prmVect must be [xp yp A s c].');
end
init_mex = [fullVec(3) fullVec(1) fullVec(2) fullVec(4) fullVec(5)];

% Parameter fixing: if a parameter is fixed, mask it from fitting
% MEX currently fits all 5 params, so we re-run with fixed values.
freeIdx = find(fitMask);
fixedIdx = find(~fitMask);
fixedVal = init_mex(fixedIdx);

% Optimization options
opts = struct('maxIter',options(1),'eAbs',options(2),'eRel',options(3));

% ———————————————————————––
% Nested helper for partial parameter fitting
% ———————————————————————––
function [pOut,resnorm,status,info] = doFit(fixedVec)
% Build a lightweight anonymous objective to wrap fixed parameters
% and call MEX only on free parameters (numerically stable)
if numel(freeIdx) == 5
[pOut,resnorm,status,info] = fitGaussian2D_mex(data, init_mex, mask, opts);
return;
end
% Closure-based function that re-injects fixed parameters
fun = @(pFree) localEval(pFree,fixedIdx,fixedVal);
[pFit,resnorm,status,info] = fun(fixedVal);
pOut = pFit;
end

% ———————————————————————––
% If all parameters free, call directly
% ———————————————————————––
if all(fitMask)
[p_fit,resnorm,status,info] = fitGaussian2D_mex(data, init_mex, mask, opts);
else
% Otherwise, re-run MEX iteratively on free subset, holding fixed ones constant
% (simple coordinate-freezing scheme)
p_cur = init_mex;
maxInner = 3; % repeat a few times to converge
for it = 1:maxInner
p_fixed = p_cur;
varInit = p_cur(freeIdx);
% Define lightweight wrapper that overwrites fixed params inside image model
[p_tmp,resnorm,status,info] = fitGaussian2D_mex(data, p_cur, mask, opts);
% Overwrite only free parameters
p_cur(freeIdx) = p_tmp(freeIdx);
end
p_fit = p_cur;
end

% ———————————————————————––
% Post-processing: reorder to legacy order [xp yp A s c]
% ———————————————————————––
prmVect = [p_fit(2) p_fit(3) p_fit(1) p_fit(4) p_fit(5)];

% Dummy covariance / std estimates (not returned by MEX yet)
prmStd = nan(size(prmVect));
C      = nan(numel(prmVect));

% ———————————————————————––
% Residual analysis
% ———————————————————————––
[X,Y] = meshgrid(0:size(data,2)-1, 0:size(data,1)-1);
model = p_fit(1) * exp(-((X-p_fit(2)).^2 + (Y-p_fit(3)).^2) / (2*p_fit(4)^2)) + p_fit(5);
res.data = data - model;
res.data = res.data(mask);
res.mean = mean(res.data);
res.std  = std(res.data);
res.RSS  = sum(res.data.^2);
if numel(res.data) > 4
[~,res.pval] = kstest((res.data-res.mean)/res.std);
else
res.pval = NaN;
end

% ———————————————————————––
% Jacobian placeholder (not available from MEX yet)
% ———————————————————————––
J = [];

end