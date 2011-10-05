%c = computeB3SplineCoefficientsFourier(s, varargin) computes cubic smoothing spline coefficients
% This function is based on the formalism/algorithm described in:
% Unser et al., IEEE Trans. Image Proc. 41(2), pp. 834-848, 1993
% For details, see Section IV. B, Eq. 4.3

% Francois Aguet, 2010 (Last modified: 10/05/2011)

function c = computeB3SplineCoefficientsFourier(s, varargin)

ip = inputParser;
ip.addRequired('s', @isvector);
ip.addOptional('lambda', 0, @(x) isscalar(x) && x>=0);
ip.addOptional('boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic'})));
ip.parse(s, varargin{:});
lambda = ip.Results.lambda;

N = length(s);

% mirror signal
if strcmpi(ip.Results.boundary, 'symmetric')
    s = [s s(N-1:-1:2)];
    M = 2*N-2;
else % periodic
    M = N;
end

% frequency vector
w = (0:M-1)*2*pi/M;

% Fourier transform of spatial input signal
S = fft(s);

% Smoothing spline pre-filter
H = 3 ./ (2+cos(w)+6*lambda*(cos(2*w)-4*cos(w)+3));

% Spline coefficients
c = real(ifft(S.*H));
c = c(1:N);
