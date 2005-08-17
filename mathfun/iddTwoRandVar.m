function S = iddTwoRandVar(s1,s2,r1,r2,s)
%iddTwoRandVar: This function identifies the dependency between two random 
%               variables by fitting a statistic model to them.
%
% SYNOPSIS: S = iddTwoRandVar(s1,s2)
%
% INPUT:
%    s1 : A vector that stores the sampling data of the first random
%         variables.
%    s2 : A vector that stores the sampling data of the second random
%         variables.
%
% OUTPUT:
%    S : A structure that contains the identified dependency parameters 
%        in the model.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First, separate the sampling data into two subsets divided by the  %%% 
%%% median of the first random variables.                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
medS1 = median(s1);
stdS1 = std(s1);
medS2 = median(s2);
dataSet1 = find(s1 <= medS1); % & s2 <= medS2);
dataSet2 = find(s1 >= medS1); % & s2 >= medS2);
%dataSet1 = find(s1 <= medS1+stdS1/2); % & s2 <= medS2);
%dataSet2 = find(s1 <= medS1-stdS1/2); % & s2 >= medS2);
%dataSet1 = find(s2 <= medS2);
%dataSet2 = find(s2 >= medS2);
s11 = s1(dataSet1);
s12 = s1(dataSet2);
s21 = s2(dataSet1);
s22 = s2(dataSet2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assemble Data.                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D11  = mean(s11);
D12  = mean(s12);
D1   = mean(s1);
D21  = mean(s21);
D22  = mean(s22);
D2   = mean(s2);
D31  = var(s11);
D32  = var(s12);
D3   = var(s1);
D41  = var(s21);
D42  = var(s22);
D4   = var(s2);
D51  = mean((s11-D11).*(s21-D21));
D52  = mean((s12-D12).*(s22-D22));
D5   = mean((s1-D1).*(s2-D2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Some combination of data.                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D3_ = D31+D32-2*D3;
D4_ = D41+D42-2*D4;
D5_ = D51+D52-2*D5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for the unkowns in the statistic model.                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = (D22-D21)/(D12-D11);
beta  = D21-alpha*D11;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Assemble the structure for output.                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.alpha = alpha;
S.beta  = beta;
