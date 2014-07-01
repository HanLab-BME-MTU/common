function [muCC,muCI,lags,xCorr] = getAverageCorrelation(TS1,varargin)
%This function estimates the average autocorrelation or cross-correlation given TS1 or TS1,TS2 inputs matrices, respectively.
%USAGE:
%
%       for ACF: [muACF,muCI,lags,allACF]  = getAverageCorrelation(TS1,varargin)
%
%       for CCF: [muCC,muCI,lags,allxCorr] = getAverageCorrelation(TS1,'TS2',TS2,varargin)
%
%Input:
%       TS1         - matrix of time series (# of variables,# of Observations) 
%                     format based on the output of the windowing package        
%                     The number of variables is the number of windows   
%                     IMPORTANT: This is the reference signal      
%
%       TS2         - same formart as TS1. Use to estimate the cross-correlation between each row of TS1 and TS2(default:[])
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
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('TS1',@(x) isnumeric(x));
[nVar1,nObs1] = size(TS1);
ip.addOptional('TS2',       TS1,@(x) isnumeric(x));
ip.addParamValue('nBoot',   1e3,@isscalar);
ip.addParamValue('alpha',   .05,@isscalar); 
ip.addParamValue('corrType','Pearson', @ischar)
ip.addParamValue('local',    length(TS1)-1,@isscalar);
ip.addParamValue('maxLag',   0,@isscalar);
ip.addParamValue('robust',   false,@islogical);

ip.parse(TS1,varargin{:});
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

xCorr = nan(2*maxLag + 1,nVar1);
bound = nan(2,nVar1);

for iVar = 1:nVar1
    
    [auxCC,bound(:,iVar),lags]        = nanCrossCorrelation(TS1(iVar,:),TS2(iVar,:),'corrType',corrT,'maxLag',maxLag,'local',local,'robust',robust);
    adjust                            = maxLag - ((numel(auxCC) - 1)/2); 
    [xCorr(adjust+1:end-adjust,iVar)] = auxCC;
    
end

%% Bootstrapping correlations
if nVar1 > 1
    
    [muCC,muCI] = correlationBootstrap(xCorr,bound(1,:),nBoot,alpha);
    
else
    
    muCC = xCorr;
    muCI = bound;
    
end