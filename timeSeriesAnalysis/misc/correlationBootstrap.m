function [meanCC,CI] = correlationBootstrap(CC,CB,nBoot,alpha)

%This function gives an average auto or cross-correlation coeficient and the
%confidence interval 
%
%Synopsis: 
%         [meanCC,CI]=correlationBootstrap(CC,CI,nBoot,alpha)
%
%Input:
%         CC - auto or cross-correlation matrix
%                CC(# of lags,# of variables)   
%
%         CB - confidence bounds for H0: CC(lag,variable)=0
%                CB is a vector of the upper bound for each variable
%                CB(1,# of variables)
%           
%         nBoot - number of bootstraped samples
%
%         alpha - alpha used to generate the bootstrap confidence intervals  
%
%Output:
%
%        meanXC - average value for the auto or cross-correlation   
%        CI     - upper and lower confidence intervals   
%
% Marco Vilela, 12/2011

if nargin < 3
    nBoot    = 1000;
end

if nargin < 4
    alpha    = 0.05;
end

[nLag,~] = size(CC);
workCC   = CC;

%Elimination of values below the CB
workCC( abs(CC) < repmat( CB,nLag,1 ) ) = 0;

opt = statset('UseParallel','always');

[confI,statCC] = bootci(nBoot,{@(x) mean(atanh(x)),workCC'},'alpha',alpha,...
          'type','bca','Options',opt);

CI     = tanh( confI );
meanCC = tanh( mean( statCC ) );

% Arc tangent is a variance-stabilizing technique
%Book: Zoubir, Iskander. Bootstrap techniques for signal processing. Page
%53

