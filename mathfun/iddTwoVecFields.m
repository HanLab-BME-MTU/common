function S = iddTwoVecFields(Va,Vh)
%iddTwoVecFields: This function identifies the dependency between two vector 
%                 fields by fitting a statistic model to them.
%
% SYNOPSIS: S = iddTwoVecFields(Va,Vh)
%
% INPUT:
%    Va : An n-by-2 matrix that stores data of the first vector field where n
%         is the number of vectors.
%    Vh : An n-by-2 matrix that stores data of the second vector field.
%
% OUTPUT:
%    S : A structure that contains the identified parameters in the model.
%
% The Model:
%
%          Va = V       + sigma1*Ra
%          Vh = alpha*V + sigma2*Rh
%
% where
%    V      : The coupling vector between the two fields.
%    Ra     : Random vector (or noise) in the first field.
%    Rh     : Random vector (or noise) in the second field.
%    alpha  : Coupling coefficient. 
%    sigma1 : Standar deviation of 'Ra'.
%    sigma2 : Standar deviation of 'Rh'.
%
% Assumptions:
%    1. The mean of 'Ra' (but not 'Rh'): ERa = (0,0). 
%    2. The variance of 'Ra' and 'Rh': 
%       var(Ra) = (1,1), var(Rh) = (1,1).
%    3. The corresponding components of 'Ra' and 'Rh' are independent:
%             E(Ra1*Rh1) = 0 and E(Ra2*Rh2) = 0.
%       Note: This is equivalent to say
%             cov(Ra1,Rh1) = 0 and cov(Ra2,Rh2) = 0
%         since ERa = (0,0).
%    4. The corresponding components of 'V' and 'Ra' are independent:
%             E(V1*Ra1) = 0 and E(V2*Ra2) = 0.
%    5. The expected value of the inner dot product of 'V-EV' and 'Rh-ERh' 
%       is zero:
%             E((V-EV).*(Rh-ERh)) = 0
%       Note: This is different from 4.
%
% Unknowns (10) to be determined:
%    alpha, sigma1, sigma2 and
%    EV     : (=(EV1,EV2)) The mean (or expected value) of 'V'.
%    var(V) : (=(VV1,VV2)) The variance of 'V'.
%    ERh    : (=(ERh1,ERh2)) The mean of 'Rh'.
%    EV1Rh1 : The mean of 'V1*Rh1'.
%
%    Note: EV1Rh1 and EV2Rh2 are related by Assumption 5:
%          (EV1Rh1+EV2Rh2) - (EV1*ERh1+EV2*ERh2) = 0.

% Data:
%    D1  : (= EVa1) The mean of 'Va1'.
%    D2  : (= EVa2) The mean of 'Va2'.
%    D3  : (= VVa1) The variance of 'Va1'.
%    D4  : (= VVa2) The variance of 'Va2'.
%    D5  : (= EVh1) The mean of 'Vh1'.
%    D6  : (= EVh2) The mean of 'Vh2'.
%    D7  : (= VVh1) The variance of 'Vh1'.
%    D8  : (= VVh2) The variance of 'Vh2'.
%    D9  : (= EVa1Vh1) The mean of 'Va1*Vh1'.
%    D10 : (= EVa2Vh2) The mean of 'Va2*Vh2'.

%%% Assemble Data.
D1  = mean(Va(:,1));
D2  = mean(Va(:,2));
D3  = var(Va(:,1));
D4  = var(Va(:,2));
D5  = mean(Vh(:,1));
D6  = mean(Vh(:,2));
D7  = var(Vh(:,1));
D8  = var(Vh(:,2));
D9  = mean(Va(:,1).*Vh(:,1));
D10 = mean(Va(:,2).*Vh(:,2));

%%% Solve the 10 equations that relate data to unknowns.
%'alpha' is the first parameter to be solved from a quadratic equation:
%    c1*alpha^2 + c2*alpha + c3 = 0
%where
%    c1 = D4-D3
%    c2 = 2*(D9-D10-D1*D5+D2*D6)
%    c3 = D8-D7
%
c1 = D4-D3;
c2 = 2*(D9-D10-D1*D5+D2*D6);
c3 = D8-D7;
if c1 == 0
   if c2 == 0
      error(['The equation can not be solved. Something is wrong ' ...
         'with the velocity based statistic model.']);
   else
      alpha = -c3/c2;
   end
else
   alpha(1) = (-c2+sqrt(c2^2-4*c1*c3))/2/c1;
   alpha(2) = (-c2-sqrt(c2^2-4*c1*c3))/2/c1;
end

%Then, we get 'VV1' and 'VV2'.
% VV1 = (alpha*(D3-D4)+D10-D2*D6-D1*D5+D9)./alpha/2;
% VV2 = (D9+D10-(D1*D5+D2*D6)-alpha*(D3-D4))./alpha/2;

%Then, we get 'VV1+VV2'.
VV1plus2 = (D10+D9-D2*D6-D1*D5)./alpha;

%Next, 'sigma1' and 'sigma2'. We calculate the square of 'sigma1' and 
% 'sigma2' which is the variance first.
% sigma1p2 = D3-VV1;
% sigma2p2 = D7-alpha.^2.*VV1-2*alpha.*(D9-alpha.*(VV1+D1^2)- ...
%    D1*(D5-alpha*D1));

sigma1p2 = ((D3+D4) - VV1plus2)/2;
sigma2p2 = ((D7+D8) - alpha.^2.*VV1plus2)/2;

j = 1;
for k = 1:length(sigma1p2)
   if sigma1p2(k) >= 0 & sigma2p2(k) >= 0
      validSol(j) = k;
      j=j+1;
   end
end

%Now, we calculate 'VV1' and 'VV2'
VV1 = D3 - sigma1p2(validSol);
VV2 = D4 - sigma1p2(validSol);

if isempty(validSol)
   error(['There is no valid solution. Something is wrong either ' ...
      'with solving the system or with the model.']);
end

sigma1 = sqrt(sigma1p2(validSol));
sigma2 = sqrt(sigma2p2(validSol));

%The rest unknowns are easy.
EV1 = D1;
EV2 = D2;

ERh1 = (D5-alpha(validSol)*D1);
ERh2 = (D6-alpha(validSol)*D2);

EV1Rh1 = (D9-alpha(validSol).*(VV1(validSol)+D1^2));
EV2Rh2 = (EV1.*ERh1+EV2.*ERh2) - EV1Rh1;

%Assemble the identified parameters into the output structure.
S.alpha  = alpha(validSol);
S.sigma1 = sigma1;
S.sigma2 = sigma2;
S.EV1    = EV1;
S.EV2    = EV2;
S.VV1    = VV1(validSol);
S.VV2    = VV2(validSol);
S.ERh1   = ERh1;
S.ERh2   = ERh2;
S.EV1Rh1 = EV1Rh1;
S.EV2Rh2 = EV2Rh2;
