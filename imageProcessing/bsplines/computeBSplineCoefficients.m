%c = computeBSplineCoefficients(s, varargin)

% Francois Aguet, 06/21/2011

function c = computeBSplineCoefficients(s, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('s', @isnumeric);
ip.addOptional('n', 3, @(x) ismember(x, [1 2 3]));
ip.addOptional('boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic', 'replicate', 'zeros'})));
ip.parse(s, varargin{:});


switch ip.Results.n
    case 3
        z1 = -2+sqrt(3);
        c0 = 6;
    case 2
        z1 = -3+2*sqrt(2);
        c0 = 8;
end

if ismember(ip.Results.n, [2 3])

    N = length(s);
    cp = zeros(1,N);
    cn = zeros(1,N);
    c = zeros(1,N);
    
    % Recursively compute coefficients
    cp(1) = getCausalInitValue(s, z1, ip.Results.boundary);
    for k=(2:N)
        cp(k) = s(k) + z1*cp(k-1);
    end
    cn(N) = getAntiCausalInitValue(cp, z1, ip.Results.boundary);
    for k=(N-1:-1:1)
        cn(k) = z1*(cn(k+1)-cp(k));
    end
    for k=(1:N)
        c(k) = c0*cn(k);
    end
else % n = 1
    c = s;
end



function c0 = getCausalInitValue(s, a, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('s', @isnumeric);
ip.addRequired('a', @isscalar);
ip.addOptional('boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic', 'replicate', 'zeros'})));
ip.parse(s, a, varargin{:});

switch ip.Results.boundary
    case 'symmetric'
        N = size(s,2);
        k = repmat(0:2*N-3, [size(s,1) 1]);
        s = [s s(:,end-1:-1:2)];
        c0 = sum(s.*a.^k, 2) / (1 - a^(2*N-2));
    case 'periodic'
        N = size(s,2);
        k = repmat(0:N-1, [size(s,1) 1]);
        s = [s(:,1) s(:,N:-1:2)];
        c0 = sum(s.*a.^k)/(1-a^N);
    case 'zeros'
        c0 = s(:,1);
    case 'replicate'
        c0 = s(:,1)/(1-a);
end



function c0 = getAntiCausalInitValue(c, a, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('c', @isnumeric);
ip.addRequired('a', @isscalar);
ip.addOptional('boundary', 'symmetric', @(x) any(strcmpi(x, {'symmetric', 'periodic', 'replicate', 'zeros'})));
ip.parse(c, a, varargin{:});

switch ip.Results.boundary
    case 'symmetric'
        N = size(c,2);
        c0 = (a/(a*a-1))*(c(:,N)+a*c(:,N-1));
    case 'periodic'
        N = size(c,2);
        %c0 = (a/(a*a-1))*(c(N)+a*c(1));
        c = [c(:,N) c c(:,1:N-3)];
        k = repmat(0:2*N-3, [size(c,1) 1]);
        c0 = -a/(1-a^N)*sum(a.^k.*c,2);
    case 'zeros'
        c0 = -a*c(:,end);
    case 'replicate'
        c0 = -c(:,end)*a/(1-a);
end

% function c0 = getAntiCausalInit_sum(s, a)
% N = length(s);
% s = [s(N) s(N-1:-1:1) s(2:N-1)];
% k = (0:2*N-3);
% c0 = -a/(1-a^(2*N-2))*sum(a.^k.*s);
