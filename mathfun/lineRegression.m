function [rhoAlpha,slopeIntcpt,dtls]=lineRegression(data,sData,rhoAlpha0,sApprox,sigma)
%LINEREGRESSION computes the model parameters of a line according to total LS
%
% SYNOPSIS [rhoAlpha,slopeIntcpt,dtls]=lineRegression(data,sData,rhoAlpha0,sApprox,sigma)
%
% INPUT data : 2xn matrix with the point coordinates
%       sData: either [] if no uncertainty is given
%                     [sData_1, sData_2] to define one global uncertainty 
%                     2xn matrix with an uncertainty value for all data points
%       rhoAlpha0: either [] if no approxation is given 
%                  [rho0, alpha0]
%       sApprox :  uncertainty of approximations if given;
%                  if rhoAlpha0 ~= [] but sApprox == [], then sApprox is set 
%                  to a large value.
%       sigma :    stdev of unit weight obervation
%
% OUTPUT  rhoAlpha : estimate of the line parameters
%         slopeIntcpt : parameters tranform3ed to slope and intercept with x_1 = 0;
%         dtls : data structure with the field
%            *.sigma  : aposteriori stdev of unit weight observation
%            *.sRhoAlpha : precision of the parameters propagated through TLS
%            *.resid  : 2xn matrix with coordinate residuals

nPts = size(data,2);

% check the available information
if(isempty(sData))
   sData = ones([2,nPts])*sigma;
else
   if(size(sData,1)*size(sData,2) == 2)
      aux(1,:) = ones([1,nPts]).*sData(1);
      aux(2,:) = ones([1,nPts]).*sData(2);
      sData = aux;
   end;
end;

if(isempty(rhoAlpha0))
   rhoAlpha0 = getApprox(data);
   sApprox = [10^32,10^32];
else
   if(isempty(sApprox))
      sApprox = [10^32,10^32];
   end;
end;

[dtls.resid,rhoAlpha,dtls.sRhoAlpha,dtls.sigma,dtls.nIter] = ...
   mexLineRegression(data',sData',rhoAlpha0,sApprox,sigma);

if(rhoAlpha(2) ~= 0)
   slopeIntcpt(1)=-cos(rhoAlpha(2))/sin(rhoAlpha(2));
   slopeIntcpt(2)=rhoAlpha(1)/sin(rhoAlpha(2));
else
   slopeIntcpt = [];
end;

dtls.resid = dtls.resid';

%--------------------------------------------------------------------------
function [ra] = getApprox(data)
% determines the parameters over the eigenvalues of the second momentum

m = mean(data,2);
mat = zeros(2);
for(i = 1:size(data,2))
   mat = mat + (data(:,i)-m) * (data(:,i) - m)';
end;

[V,D] = eig(mat);

if(D(1,1) < D(2,2))
   ra(2) = mod(atan2(V(2,1),V(1,1))+pi,pi);
else
   ra(2) = mod(atan2(V(2,2),V(1,2))+pi,pi);
end;
ra(1) = m' * [cos(ra(2));sin(ra(2))];

    
   