% ADTEST(X) implements the Anderson-Darling normality test, assuming unknown mean and variance
%
% Inputs:
%          x : Sample vector (test only work for n>=5)
%    {alpha} : Significane level. Possible values: 0.01 0.025 0.05 (default) 0.1 0.15
%
% Output:
%          H : Result of the hypothesis test
%              0: Do not reject H0 at given significance level
%              1: Reject H0 at given significance level
%         A2 : Adjusted test statistic
%       cval : Critical value for the test
%
%
% Critical values taken from 
% [1] M.A. Stephens, J. Am. Stat. Assoc. 69(347), pp. 730-737, 1974
%
% See also
% [2] http://en.wikipedia.org/wiki/Anderson-Darling_test
%
% Note: this function only implements the case where both the mean and variance
% of the sample are unknown.

% Francois Aguet, 01/28/2012

function [H, A2, cval] = adtest(x, varargin)

alphaVec = [0.01 0.025 0.05 0.1 0.15];

ip = inputParser;
ip.addRequired('x');
ip.addOptional('alpha', 0.05, @(x) ismember(x, alphaVec));
ip.addParamValue('mu', []);
ip.addParamValue('sigma', []);
ip.parse(x, varargin{:});

n = numel(x);
if n<5
    error('The Anderson-Darling test requires at least 5 samples.');
end
x = reshape(x, [1 n]);

estidx = [0 0]; % [sigma mu]

mu = ip.Results.mu;
if isempty(mu)
    mu = mean(x);
    estidx(2) = 1;
end

sigma = ip.Results.sigma;
if isempty(sigma)
    sigma = std(x);
    estidx(1) = 1;
end

% case
c = bin2dec(num2str(estidx))+1;

% Case 1: mean and variance known
% Case 2: mean unknown, variance known
% Case 3: mean known, variance unknown
% Case 4: mean and variance unknown

% Look-up table for critical values
ctable = [3.857 3.070 2.492 1.933 1.610;
          1.573 1.304 1.105 0.908 NaN;
          3.690 2.904 2.323 1.760 NaN;
          1.092 0.918 0.787 0.656 0.576];

      
% Alternative lookup table (d'Agostino 1986)
% ctable(4,:) = [1.035 0.873 0.752 0.631 NaN];
      
cval = ctable(c,alphaVec==ip.Results.alpha);
      
% Test statistic
x = sort(x);
z = normcdf(x, mu, sigma);
i = 1:n;
A2 = -n - 1/n*sum( (2*i-1).*(log(z) + log(1-z(n+1-i))) );

% correction factor for sample mean and variance
if c==4
    A2 = A2 * (1+4/n-25/n^2);
    % A2 = A2 * (1.0 + 0.75/n + 2.25/n^2); % for the table by d'Agostino
end
H = A2 > cval;
