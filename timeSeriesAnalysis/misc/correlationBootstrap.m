function [meanCC,CI] = correlationBootstrap(CC,CB,varargin)
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
%Optional input:
%
%         nBoot - number of bootstraped samples (default value 1000)
%
%         alpha - alpha used to generate the bootstrap confidence intervals
%         (default value 0.05)
%
%Output:
%
%        meanXC - average value for the auto or cross-correlation   
%        CI     - upper and lower confidence intervals   
%
% Marco Vilela, 12/2011

%Input check
%SB: still need to find minimimum values for number of samples/variables
ip=inputParser;
ip.addRequired('CC',@(x) isnumeric(x) && size(x,2)>2);
[nLag,nVar] = size(CC);
ip.addRequired('CB',@(x) isnumeric(x) && isequal(size(x),[1 nVar]));
ip.addOptional('nBoot',1e3,@isscalar);
ip.addOptional('alpha',.05,@isscalar);
ip.parse(CC,CB,varargin{:});
nBoot=ip.Results.nBoot;
alpha=ip.Results.alpha;

% Fix bug when correlation function is equal to 1 or -1
workCC   = CC;
workCC(abs(CC)==1)=workCC(abs(CC)==1)-sign(workCC(abs(CC)==1))*1e-9;

%Elimination of values below the CB
workCC( abs(CC) < repmat( CB,nLag,1 ) ) = 0;

opt = statset('UseParallel','always');

[confI,statCC] = bootci(nBoot,{@(x) nanmean(atanh(x)),workCC'},'alpha',alpha,...
          'type','bca','Options',opt);

CI     = tanh( confI );
meanCC = tanh( mean( statCC ) );

% Arc tangent is a variance-stabilizing technique
%Book: Zoubir, Iskander. Bootstrap techniques for signal processing. Page
%53

