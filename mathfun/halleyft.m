function [ refined ] = halleyft( v, guess, freq, deriv, TOL, maxIter, varargin )
%newtonft Does Newton iteration to refine roots of Fourier series
%
% INPUT
% v - known values of Fourier series at regular intervals or it's Fourier
% transform. Each Fourier series is encoded in the first dimension as a
% column
%
% guess - guess for the root to find, see interpft_extrema. Multiple
% guesses for the same Fourier series are encoded in the first dimension as
% a column. All other dimensions must match v.
%
% freq - logical. If true, then v is the Fourier transform of the values
%
% deriv - solve for zero of the derivative indicated
%         (optional, default = 0);
%
% TOL - tolerance
%
% OUTPUT
% refined - refined zero

if(nargin < 3)
    freq = false;
elseif(ischar(freq))
    ip = inputParser;
    ip.addParameter('freq',false);
    ip.addParameter('deriv',0);
    ip.addParameter('TOL',1e-12);
    ip.addParameter('maxIter',10);
    argsIn = {freq};
    if(nargin > 3)
        argsIn = [argsIn deriv];
    end
    if(nargin > 4)
        argsIn = [argsIn TOL];
    end
    if(nargin > 5)
        argsIn = [argsIn maxIter];
    end
    argsIn = [argsIn varargin];
    ip.parse(argsIn{:});
    freq = ip.Results.freq;
    deriv = ip.Results.deriv;
    TOL = ip.Results.TOL;
    maxIter = ip.Results.maxIter;
else
    if(nargin < 4)
        deriv = 0;
    end
    if(nargin < 5)
        TOL = 1e-12;
    end
    if(nargin < 6)
        % If more than 10, probably should use interpft_extrema
        maxIter = 10;
    end
end

derivs = [0 1 2] + deriv;

    %% From interpft1_derivatives
    K = floor(size(v,1)/2);
    freqM = ifftshift(-K:K).'*1i;
    if(~freq)
        v_hat = fft(v);
    else
        v_hat = v;
    end
    derivDim = ndims(v)+1;
    freqMs = arrayfun(@(x) freqM.^x,derivs,'UniformOutput',false);
    freqMs = cat(derivDim,freqMs{:});
    
    v_hat = bsxfun(@times,v_hat,freqMs);  
    
    xqrep = ones(1,derivDim);
    xqrep(derivDim) = length(derivs);
    
    colons = {':'};
    colons = colons(ones(derivDim,1));
    zeroth_d = colons;
    zeroth_d{derivDim} = 1;
    first_d = colons;
    first_d{derivDim} = 2;
    second_d = colons;
    second_d{derivDim} = 3;

    %% Perform Newton iteration
numIter = 0;
% do while
while(~numIter || any(exceeds_tol) && any(new_guess_is_better(:)))
%     disp('hi');
    guess_vals = real(interpft1([0 2*pi],v_hat,repmat(guess,xqrep),'horner_freq'));
    new_guess = guess - 2*guess_vals(zeroth_d{:}).*guess_vals(first_d{:})./(2*guess_vals(first_d{:}).^2-guess_vals(zeroth_d{:}).*guess_vals(second_d{:}));
    new_guess = wraparoundN(new_guess,0,2*pi);
    new_guess_vals = interpft1([0 2*pi],v_hat,repmat(new_guess,xqrep),'horner_freq');
    new_guess_is_better = abs(new_guess_vals(zeroth_d{:})) < abs(guess_vals(zeroth_d{:}));
    guess(new_guess_is_better) = new_guess(new_guess_is_better);
    % Establish new derivative value
    new_guess_vals = new_guess_vals(zeroth_d{:});
    guess_vals = guess_vals(zeroth_d{:});
    zero_vals = min(abs(guess_vals),abs(new_guess_vals));
    exceeds_tol = zero_vals(:) > TOL;
    numIter = numIter + 1;
    if(numIter > maxIter)
%         warning('halleyft:maxIter','halleyft: maximum iteration reached');
        break;
    end
end

guess(exceeds_tol) = NaN;

refined = guess;

end
