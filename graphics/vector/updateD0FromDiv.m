function d0=updateD0FromDiv(divM,d0in,alpha,nPi,nPg)
% updateD0FromDiv creates the matrix d0 to be used with the interpolator
%
% SYNOPSIS   d0=updateD0FromDiv(divM,d0,alpha,nPi,nPg)
%
% INPUT      divM : divergence as calculated by vectorFieldDiv.
%            d0in : initial d0. It can be a scalar or a (nPg x nPi) matrix.
%            alpha: see OUTPUT. alpha SHOULD be 0 < alpha <= 1, but may also 
%                   be set higher.
%            nPi  : number of vectors prior of interpolation (raw data)
%                   (see vectorFieldInterp).
%            nPg  : number of grid points for interpolation
%                   (see vectorFieldInterp).
%
% OUTPUT    d0    : parameter d0 for the interpolator 
%                   (see vectorFieldInterpolator).
%                   At each position, d0 will be equal to 
%                   d0/(1+alpha/max(divM(all))*divM),
%                   where alpha is a parameter defining the range of
%                   variation possible for divM.
%                   d0 has size (nPg x nPi).

% Check input parameters
if nargin~=5
    error('The function expects 5 input parameters.');
end
if size(divM,1)~=nPg
    error('Parameter nPg does not match the number of rows of divM');
end
if size(d0in)~=[1 1] & size(d0in)~=[nPg nPi]
    error('The dimensions of d0 are wrong');
end
if alpha<=0
    error('The parameter alpha must be positive');
end

% Dimensions of the input d0
[d0y d0x]=size(d0in);

% Initialize d0
d0=zeros(nPg,nPi);

% Fill d0 to be used in the correlation matrix for the interpolator
if d0y==1
    for i=1:nPg
        d0(i,1:nPi)=d0in/(1+alpha/max(divM(:,3))*abs(divM(i,3)));
    end
else
    for i=1:nPg
        d0(i,1:nPi)=d0in(i,1:nPi)/(1+alpha/max(divM(:,3))*abs(divM(i,3)));
    end
end
