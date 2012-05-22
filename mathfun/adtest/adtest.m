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
% For the test and its derivation, see
% [1] Anderson & Darling, Ann. Math. Stat. 23, 1952
% [2] Anderson & Darling, J. Am. Stat. Assoc. 49, 1954
%
% Critical values taken from 
% [3] M.A. Stephens, J. Am. Stat. Assoc. 69(347), pp. 730-737, 1974
% [4] M.A. Stephens, Ann. Stat. 4(2), pp. 357-369, 1976
%
% See also
% [5] http://en.wikipedia.org/wiki/Anderson-Darling_test

% Francois Aguet (last modified 03/27/2012)

function [H, pval, A2, cval] = adtest(x, varargin)

alphaVec = [0.01 0.025 0.05 0.1 0.15];

ip = inputParser;
ip.addRequired('x');
ip.addOptional('alpha', 0.05, @(x) ismember(x, alphaVec));
ip.addParamValue('mu', []);
ip.addParamValue('sigma', []);
ip.addParamValue('Distribution', 'normal', @(x) any(strcmpi(x, {'normal', 'exponential'})));
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


% sort samples in ascending order
x = sort(x);

% case
switch ip.Results.Distribution
    case 'normal'
        c = bin2dec(num2str(estidx))+1;
        z = normcdf(x, mu, sigma);

    case 'exponential'
        c = 5;
        z = expcdf(x, mu);
end

% Look-up table for critical values
%  Case 1: mean and variance known
%  Case 2: mean unknown, variance known
%  Case 3: mean known, variance unknown
%  Case 4: mean and variance unknown
%  Case 5: exponential distribution, mean unknown

% Values from [3]
ctable = [3.857 3.070 2.492 1.933 1.610;
          1.573 1.304 1.105 0.908 NaN;
          3.690 2.904 2.323 1.760 NaN;
          1.092 0.918 0.787 0.656 0.576
          1.957 1.606 1.341 1.078 0.922];

% Values from [4]
% ctable = [3.857 3.070 2.492 1.933 1.610;
%           1.541 1.281 1.088 0.897 0.784;
%           3.682 2.890 2.315 1.761 1.443;
%           1.029 0.870 0.751 0.632 0.560;
%           1.943 1.587 1.326 1.070 0.918];

% Another alternative lookup table (d'Agostino 1986)
% ctable(4,:) = [1.035 0.873 0.752 0.631 NaN];


cval = ctable(c,alphaVec==ip.Results.alpha);
      
% Test statistic
i = 1:n;
A2 = -n - 1/n*sum( (2*i-1).*(log(z) + log(1-z(n+1-i))) );

% correction factor for sample mean and variance
pval = NaN;
if c==4
    %A2 = A2 * (1+4/n-25/n^2);
    A2 = A2 * (1.0 + 0.75/n + 2.25/n^2); % for the table by d'Agostino
    if (0.600<A2 && A2<=13)
        pval = exp(1.2937 - 5.709*A2 + 0.0186*A2*A2);
    end
    if (0.340<A2 && A2<=0.600)
        pval = exp(0.9177 - 4.279*A2 - 1.38*A2*A2);
    end
    if (0.200<A2 && A2<=0.340)
        pval = 1 - exp(-8.318 + 42.796*A2 - 59.938*A2*A2);
    end
    if (A2<=0.200)
        pval = 1 - exp(-13.436 + 101.14*A2 - 223.73*A2*A2);
    end
end
H = A2 > cval;
