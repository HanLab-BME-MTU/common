function [ Y ] = idct1( X, n, dim, enforceSymmetry )
%idct1 Inverse discrete cosine transform Type I
%
% Same as dct1 except it multiplies the result by 2/(n-1)
%
% See dct1, idct
%
% Mark Kittisopikul, January 26th, 2017
% Lab of Khuloud Jaqaman
% UT Southwestern

% Default to dim = 1.
if(nargin < 3 || isempty(dim))
    dim = 1;
    % If row vector, then dim = 2
    if(isrow(X))
        dim = 2;
    end
end

if(nargin < 2 || isempty(n))
    n = size(X,dim);
end

if(nargin < 4)
    enforceSymmetry = false;
end

% Normalize by 2/(n-1) to invert the dct1
Y = dct1(X,n,dim,enforceSymmetry)*2/(n-1);


end

