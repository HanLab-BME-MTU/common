function [ refined, refined_derivs ] = halleyft( v, guess, freq, deriv, TOL, maxIter, avoidNaN )
%newtonft Does Halley iteration to refine roots of Fourier series
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
% TOL - tolerance for absolute distance from 0, default: 1e-12
%
% maxIter - maximum number of iterations, default: 10
%
% avoidNaN - For non-vector input, avoid processing NaN initial guesses.
%            The matrix should be sorted for this to work.
%
% OUTPUT
% refined - refined zero

if(nargin < 3)
    freq = false;
end
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
if(nargin < 7)
    avoidNaN = isvector(guess);
end

derivs = [0 1 2] + deriv;

    %% From interpft1_derivatives
    if(~freq)
        v_hat = fft(v(:,:));
    else
        v_hat = v(:,:);
    end
    
    if(avoidNaN)
        % NaN guesses slow down the algorithm
        % Group columns by number of non-nan guesses and recurse
        % Depends on the guess being sorted such that NaNs are all in the
        % bottom iterations
        guess_sz = size(guess);
        guess = guess(:,:);
        guess_count = size(guess,1)-sum(isnan(guess),1);
        refined = guess;
        if(nargout > 1)
            refined_derivs = guess(:,:,[1 1 1]);
            for i=1:size(guess,1)
                s = guess_count == i;
                [refined(1:i,s),refined_derivs(1:i,s,:)] = halleyft(v_hat(:,s), guess(1:i,s), true, deriv, TOL, maxIter, false);
            end
        else
            for i=1:size(guess,1)
                s = guess_count == i;
                refined(1:i,s) = halleyft(v_hat(:,s), guess(1:i,s), true, deriv, TOL, maxIter, false);
            end
        end
        refined = reshape(refined,guess_sz);
        return;
    end
    
    K = floor(size(v,1)/2);
    freqM = ifftshift(-K:K).'*1i;
    
%     derivDim = ndims(v)+1;
    derivDim = 3;
    freqMs = arrayfun(@(x) freqM.^x,derivs,'UniformOutput',false);
    freqMs = cat(derivDim,freqMs{:});
    
    v_hat = bsxfun(@times,v_hat,freqMs);  
    
    xqrep = ones(1,derivDim);
    xqrep(derivDim) = length(derivs);
    

    %% Perform iteration
args = {};
if(isa(guess,'gpuArray'))
    args = { 'gpuArray'};
end

guess_sz = size(guess);
guess = guess(:,:);
if(nargout > 1)
    refined_derivs = NaN([size(guess) 3]);
end

% Indicates that a column still needs further iterations
columnNotDone = true(1,size(guess,2),args{:});
% Tracks whether a new guess moves the derivative closer to zero
new_guess_is_better = true(size(guess),args{:});

% Calculate initial value
guess_vals = interpft1([0 2*pi],v_hat,repmat(guess,xqrep),'horner_freq');
% do while
% while(~numIter || any(exceeds_tol) && any(new_guess_is_better(:)))
numIter = 0;
while(~numIter || any(columnNotDone))
    disp('hi');
    v = guess_vals(:,:,1);
    vd = guess_vals(:,:,2);
    vdd = guess_vals(:,:,3);
    new_guess = guess(:,columnNotDone) - 2*v.*vd./(2*vd.^2-v.*vdd);
    new_guess = wraparoundN(new_guess,0,2*pi);
    new_guess_vals = interpft1([0 2*pi],v_hat(:,columnNotDone,:),repmat(new_guess,xqrep),'horner_freq');
    new_guess_is_better_now = abs(new_guess_vals(:,:,1)) < abs(guess_vals(:,:,1));
    new_guess_is_better(:,columnNotDone) = new_guess_is_better_now;
    guess(new_guess_is_better) = new_guess(new_guess_is_better_now);
    if(nargout > 1)
        refined_derivs(new_guess_is_better(:,:,[1 1 1])) = new_guess_vals(new_guess_is_better_now(:,:,[1 1 1]));
    end
%     columnNotDone(columnNotDone) = any(new_guess_is_better(:,columnNotDone),1);

    % Establish new derivative value
    zero_vals = min(abs(guess_vals(:,:,1)),abs(new_guess_vals(:,:,1)));
    exceeds_tol = zero_vals > TOL;
       
    numIter = numIter + 1;
    if(numIter > maxIter)
        % Maximum number of iterations exceeded
        temp = guess(:,columnNotDone);
        % Set values whose derivative exceed tolerance to NaN
        temp( exceeds_tol ) = NaN;
        guess(:,columnNotDone) = temp;
        if(nargout > 1)
            temp = refined_derivs(:,columnNotDone,:);
            temp( exceeds_tol(:,:,[1 1 1]) ) = NaN;
            refined_derivs(:,columnNotDone,:) = temp;
        end
        % Issue a warning, use warning('off','halleyft:maxIter') to
        % deactivate
        warning('halleyft:maxIter','halleyft: maximum iteration reached');
        break;
    else
        % If new guess is not better and the value still exceeds tolerance,
        % then set the guess to NaN
        temp = guess(:,columnNotDone);
        invalid = ~new_guess_is_better_now & exceeds_tol;
        temp( invalid ) = NaN;
        guess(:,columnNotDone) = temp;
        if(nargout > 1)
            temp = refined_derivs(:,columnNotDone,:);
            temp( invalid(:,:,[1 1 1]) ) = NaN;
            refined_derivs(:,columnNotDone,:) = temp;
        end

        % New guess will only be considered further if derivative exceeds
        % tolerance
        new_guess_is_better_now = new_guess_is_better_now & exceeds_tol;
        new_guess_is_better(:,columnNotDone) = new_guess_is_better_now;
        % Find columns that are not done in the set of not done columns
        newNotDone = any(new_guess_is_better_now,1);
        % Use new guess values in the next iteration
        guess_vals = new_guess_vals(:,newNotDone,:);
        columnNotDone(columnNotDone) = newNotDone;
    end
end

refined = reshape(guess,guess_sz);
if(nargout > 1)
    refined_derivs = reshape(refined_derivs,[guess_sz 3]);
end

end
