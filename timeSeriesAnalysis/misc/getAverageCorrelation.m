function [muCC,muCI,lags,xCorr,excluded] = getAverageCorrelation(TS1,varargin)
%This function estimates the average autocorrelation or cross-correlation given TS1 or TS1,TS2 inputs matrices, respectively.
%USAGE:
%
%       for ACF:        [muACF,muCI,lags,allACF] = getAverageCorrelation(TS1,varargin)
%
%       for CCF: [muCC,muCI,lags,xCorr,excluded] = getAverageCorrelation(TS1,'TS2',TS2,varargin)
%
%Input:
%       TS1 - matrix of time series (# of variables,# of Observations) 
%            format based on the output of the windowing package        
%            The number of variables is the number of windows   
%
%       TS2         - same formart as TS1. Use to estimate the cross-correlation between each row of TS1 and TS2(default:[])
%       excludeWin1 - variables(or windows) from TS1 to be excluded from the calculations(default:[])
%       excludeWin2 - same as excludeWin1 but for TS2(default:[])
%       nBoot       - number of bootstrap samples(default:1e3)
%       alpha       - alpha value used to estimate the confidence interval for the average correlation(default:0.05)
%       corrTyep    - correlation type: Pearson, Spearman or Kendall
%       local       - if corrType = Kendall, local defines a local neighborhood where the correlation is calculated
%       maxLag      - scalar. maximum correlation lag(default:0)
%       robust      - logical. Cross-correlation is calculated using robust regression if this parameter is true.(default:false)
%
%Output:
%       muCC  - Average ACF or CCF
%       muCI  - confidence interval for muCC given the input alpha
%       lags  - lags for muCC
%       xCorr - ACF or CCF for each row of TS1/TS2 
%       excluded - variables(windows) excluded = union(excludeWin1,excludeWin2)
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('TS1',@(x) isnumeric(x));
[nVar1,nObs1] = size(TS1);

ip.addOptional('excludeWin1',     [],@(x) isnumeric(x));
ip.addOptional('TS2',     TS1,@(x) isnumeric(x));
ip.addOptional('excludeWin2',     [],@(x) isnumeric(x));

ip.addParamValue('nBoot',   1e3,@isscalar);
ip.addParamValue('alpha',   .05,@isscalar);
ip.addParamValue('corrType','Pearson', @ischar)
ip.addParamValue('local',length(TS1)-1,@isscalar);
ip.addParamValue('maxLag',0,@isscalar);
ip.addParamValue('robust',false,@islogical);

ip.parse(TS1,varargin{:});
exclude1 = ip.Results.excludeWin1;
exclude2 = ip.Results.excludeWin2;
nBoot    = ip.Results.nBoot;
alpha    = ip.Results.alpha;
TS2      = ip.Results.TS2;
corrT    = ip.Results.corrType;
local    = ip.Results.local;
maxLag   = ip.Results.maxLag;
robust   = ip.Results.robust;


%% Setting TS2 up

[nVar2,nObs2] = size(TS2);
if nVar1 ~= nVar2
    error('Different number of windows(variables)')
end

if nObs1 ~= nObs2
    error('Different number of time points')
end

%% Calculating pair-wise correlations
excluded        = union(exclude2,exclude1);
TS1(excluded,:) = [];
TS2(excluded,:) = [];

nVar  = size(TS1,1);
xCorr = nan(2*maxLag + 1,nVar);
bound = nan(2,nVar);

for iVar = 1:nVar
    
    [xCorr(:,iVar),bound(:,iVar),lags] = nanCrossCorrelation(TS1(iVar,:),TS2(iVar,:),'corrType',corrT,'maxLag',maxLag,'local',local,'robust',robust);
    
end

%% Bootstrapping correlations
if nVar1 > 1
    
    [muCC,muCI] = correlationBootstrap(xCorr,bound(1,:),nBoot,alpha);
    
else
    
    muCC = xCorr;
    muCI = bound;
    
end