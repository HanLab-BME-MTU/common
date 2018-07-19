function [ Y ] = dct1( X, n, dim, enforceSymmetry )
%DCT1 Perform a Type I discrete cosine transform
%
% X is the signal of interest
% n is the degree of the dct to do, will truncate or pad with zeros
% dim - dimension to perform the transform (default: 1)
%       If a row vector is given, then dim = 2
% enforceSymmetry - ensures that the output is exactly even or odd if the
% has all odd or even values as zero
%
% This maps the problem to a discrete fourier transform. For a sequence
% abcde the dct1 is the equivalent of fft(abcdedcb).
%
% This differs from MATLAB builtin function dct which performs a Type II
% dct
%
% See also dct, chebfun.dct, fft
%
% Mark Kittisopikul, January 26th, 2017
% Lab of Khuloud Jaqaman
% UT Southwestern

colons = {':'};
colons = colons(ones(1,ndims(X)-1));

% Default to dim = 1.
if(nargin < 3 || isempty(dim))
    dim = 1;
    unshift = 0;
    % If row vector, then dim = 2
    if(isrow(X))
        dim = 2;
        unshift = -1;
    end
end

% If dim ~= 1, then shift so dim is effective 1
if(dim ~= 1)
    X = shiftdim(X,dim - 1);
    unshift = -dim+1;
end

% Default size is the dimension along dim
if(nargin < 2 || isempty(n))
    n = size(X,1);
elseif(n ~= size(X,1))
    % Pad with zeros to n, does nothing if n < size(X,1)
    X(end+1:n,colons{:}) = 0;
    % Truncate if n is shorter
    X = X(1:n,colons{:});
end

% Do not enforce symmetry by default
if(nargin < 4 || isempty(enforceSymmetry))
    enforceSymmetry = false;
end

if(enforceSymmetry)
    % If all the even cofficients are zero, make the output odd (below)
    isEven = all(X(2:2:end,:) == 0);
    % If all the odd coefficients are zero, make the output even (below)
    isOdd = all(X(1:2:end,:) == 0);
end

% Replicate with a mirror image, do not repeat the boundaries
X = [ X ; X(end-1:-1:2,colons{:}) ]/2;
Y = fft(X);
% Truncate to size, better to pass n?
Y = Y(1:n,colons{:});

% If intput is real, output should be real
if(isreal(X))
    Y = real(Y);
end

if(enforceSymmetry)
    szY = size(Y);
    % If all the even cofficients are zero (above), make the output odd
    Y(:,isEven) = (Y(:,isEven)+flipud(Y(:,isEven)))/2;
    % If all the odd coefficients are zero (above), make the output even
    Y(:,isOdd) = (Y(:,isOdd)-flipud(Y(:,isOdd)))/2;
    Y = reshape(Y,szY);
end

if(unshift)
    % Unshift so that dimensions are in the same order as input
    Y = shiftdim(Y,unshift);
end


end

