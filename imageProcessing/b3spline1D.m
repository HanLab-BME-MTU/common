% [c] = b3spline1D(signal, boundary)
% 
% Computes the cubic B-spline coefficients.
%
% Inputs: 
%           signal   : input signal
%           boundary : boundary conditions: 'mirror' (default) or 'periodic'

% Francois Aguet, June 2010

function c = b3spline1D(signal, boundary)

if nargin<2
    boundary = 'mirror';
end

% cubic spline parameters
c0 = 6;
a = -2+sqrt(3);

c = prefilter(reshape(signal, [1 length(signal)]), c0, a, boundary);
c = reshape(c, size(signal));


function cn = prefilter(signal, c0, z1, mode)
N = length(signal);
cp = zeros(size(signal));
cn = zeros(size(signal));

if (strcmp(mode, 'mirror')==1)
    cp(1) = getCausalInit_Mirror(signal, z1);
    for k = 2:N
        cp(k) = signal(k) + z1*cp(k-1);
    end;
    cn(N) = getAntiCausalInit_Mirror(cp, z1);
    for k = N-1:-1:1
        cn(k) = z1*(cn(k+1) - cp(k));
    end;
elseif (strcmp(mode, 'periodic')==1)
    cp(:,1) = getCausalInit_Periodic(signal, z1);
    for k = 2:N
        cp(k) = signal(k) + z1*cp(k-1);
    end;
    cn(:,N) = getAntiCausalInit_Periodic(cp, z1);
    for k = N-1:-1:1
        cn(k) = z1*(cn(k+1) - cp(k));
    end
else
    error('Mode: unknown boundary conditions');
end;
cn = c0*cn;


function c0 = getAntiCausalInit_Mirror(signal, a)
N = length(signal);
c0 = (a/(a*a-1))*(signal(N)+a*signal(N-1));


function c0 = getAntiCausalInit_Periodic(signal, a)
N = length(signal);
signal = [signal(N) signal signal(1:N-3)];
k = 0:2*N-3;
c0 = -a/(1-a^N)*sum(a.^k.*signal);


function out = getCausalInit_Mirror(signal, a)
N = length(signal);
k = 0:2*N-3;
signal = [signal signal(end-1:-1:2)];
out = sum(signal.*a.^k) / (1 - a^(2*N-2));


function out = getCausalInit_Periodic(signal, a)
N = length(signal);
k = 0:N-1;
signal = [signal(:,1) signal(:,N:-1:2)];
out = sum(signal.*a.^k)/(1-a^N);