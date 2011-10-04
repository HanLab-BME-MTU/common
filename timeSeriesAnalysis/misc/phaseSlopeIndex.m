function [out]=phaseSlopeIndex(X,Fs)
%This function computes the phase slope index between two time series
%Synopsis: [out]=phaseSlopeIndex(X,Fs) 
%
%Input
%X(:,1) and X(:,2)
%   Fs- Sampling frequency in Hz
%
%Output
%  If out > 0, X(:,1) comes first in relation to X(:,2), or X1 -> X2
%  If out < 0, X(:,2) comes first in relation to X(:,1), or X2 -> X1
%
%Reference
%G. Nolte, A. Ziehe, V.V. Nikulin, A. Schlogl, N. Kramer, T. Brismar, and K.R. Muller. Robustly
%estimating the flow direction of information in complex physical systems. Phys Rev Lett, 100:
%234101, 2008.
%Marco Vilela, last update 11/17/2010

if nargin < 2
    Fs=1;
end

%Each column is a variable time serie
[nobs,nvar]=size(X);
if nvar > nobs
    X=X';
    [nobs,~]=size(X);
end
%************************************

%Removing mean and linear trend******
X = X - repmat(mean(X),nobs,1);
X = detrend(X,'linear');
%************************************

Coh   = mscohere(X(:,1),X(:,2),[],[],[],Fs);%Coherence
CrPS  = cpsd(X(:,1),X(:,2),[],[],[],Fs);%Cross Spectrum
Teta  = angle(CrPS);%Cross-Spectrum angles
CpCoh = sqrt(Coh).*exp(-1i*Teta);%Complex coherence
out   = -imag(sum( conj(CpCoh(1:end-1)).*CpCoh(2:end) ) );
