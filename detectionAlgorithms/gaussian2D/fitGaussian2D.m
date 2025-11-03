function [prmVect, prmStd, C, res, J] = fitGaussian2D(data, prmVect, mode, options)
% Wrapper that calls the compiled MEX: fitGaussian2D.mex*
% Keeps original signature:
%   [prmVect prmStd C res J] = fitGaussian2D(data, prmVect, mode, options)
%
% This wrapper:
%   - sanitizes inputs (dtype/shape/NaNs)
%   - ensures defaults for 'mode' and 'options'
%   - clamps sigma to a tiny positive number
%   - retries with relaxed tolerances / larger maxIter if GSL reports no-convergence
%   - perturbs initial guess slightly and retries as a last resort
%
% If you ever need to bypass the wrapper and call the MEX from here,
% we use builtin('fitGaussian2D', ...) to avoid recursion.

    % -------- defaults & hygiene --------
    if nargin < 3 || isempty(mode),    mode = 'xyasc'; end
    if nargin < 4 || isempty(options), options = [200 1e-6 1e-6]; end
    validateattributes(data, {'numeric','logical'}, {'2d','nonempty'}, mfilename, 'data', 1);
    if ~isa(data,'double') || ~isreal(data), data = double(real(data)); end

    prmVect = prmVect(:).';
    if numel(prmVect) ~= 5
        error('fitGaussian2D:BadInit','prmVect must be [xp yp A s c].');
    end

    % Clamp sigma to > 0 to reduce numeric issues in the solver
    MIN_SIGMA = 1e-6;
    if ~(isfinite(prmVect(4))) || prmVect(4) < MIN_SIGMA
        prmVect(4) = max(MIN_SIGMA, abs(prmVect(4)));
    end

    % Make sure mode contains only valid letters
    mode = lower(char(mode));
    mode = mode(ismember(mode, 'xyasc'));  % keep only valid
    if isempty(mode)
        % if someone passed nonsense, default to all
        mode = 'xyasc';
    end

    % Ensure options has 3 elements
    options = options(:).';
    if numel(options) < 3
        options(3) = 1e-6;
    end
    if ~isfinite(options(1)) || options(1) < 1,  options(1) = 200;  end % maxIter
    if ~isfinite(options(2)) || options(2) <= 0, options(2) = 1e-6; end % eAbs
    if ~isfinite(options(3)) || options(3) <= 0, options(3) = 1e-6; end % eRel

    % -------- call MEX safely (attempts) --------
    % Attempt 1: as-is
    try
        [prmVect, prmStd, C, res, J] = builtin('fitGaussian2D', data, prmVect, mode, options);
        return;
    catch ME1
        % If the mex doesn't exist, surface a clear error
        if strcmp(ME1.identifier,'MATLAB:UndefinedFunction')
            error('fitGaussian2D:MEXNotFound', ...
                ['Could not find fitGaussian2D MEX on path.\n' ...
                 'Check with: which -all fitGaussian2D\n' ...
                 'Make sure fitGaussian2D.mexmaca64 is ahead of this wrapper.']);
        end

        % If it’s a hard error unrelated to convergence, rethrow
        if ~isGslNoConvergence(ME1)
            rethrow(ME1);
        end
        % otherwise, try relaxed settings below
    end

    % Attempt 2: relax tolerances & increase maxIter
    opt2 = options;
    opt2(1) = max(opt2(1), 1000);  % maxIter
    opt2(2) = max(opt2(2), 1e-5);  % eAbs
    opt2(3) = max(opt2(3), 1e-5);  % eRel
    % keep sigma positive
    prm2 = prmVect; prm2(4) = max(MIN_SIGMA, abs(prm2(4)));
    try
        [prmVect, prmStd, C, res, J] = builtin('fitGaussian2D', data, prm2, mode, opt2);
        return;
    catch ME2
        if ~isGslNoConvergence(ME2)
            rethrow(ME2);
        end
    end

    % Attempt 3: small perturbations to initial guess (helps when
    % starting at a flat/ill-conditioned point)
    rngState = rng; %#ok<RNGR>  % preserve user RNG
    cleaner = onCleanup(@() rng(rngState));
    rng('default');

    prm3 = prmVect;
    % tiny nudges
    prm3(1) = prm3(1) + 0.2;  % xp
    prm3(2) = prm3(2) + 0.2;  % yp
    prm3(3) = prm3(3) + 0.05*max(1,abs(prm3(3))); % A
    prm3(4) = max(MIN_SIGMA, abs(prm3(4)))*(1 + 0.05); % s
    % c는 변화 없이 두는 편이 안정적인 경우가 많음

    opt3 = opt2;
    try
        [prmVect, prmStd, C, res, J] = builtin('fitGaussian2D', data, prm3, mode, opt3);
        return;
    catch ME3
        % 마지막으로, 수렴 실패지만 best-so-far를 반환했을 가능성 경고
        warning('fitGaussian2D:NoConvergence', ...
            'GSL did not reach tolerance after retries: %s.\nReturning best-so-far if available.', ME3.message);
        rethrow(ME3); % 사용자가 try/catch로 다룰 수 있게 넘김
    end
end

function tf = isGslNoConvergence(ME)
    % Detect typical GSL no-convergence / tolerance messages coming out of the MEX
    msg = lower(ME.message);
    tf = contains(msg, 'did not reach requested tolerance') || ...
         contains(msg, 'cannot reach the specified tolerance') || ...
         contains(msg, 'no convergence') || ...
         strcmp(ME.identifier, 'fitGaussian2D:NoConvergence');
end