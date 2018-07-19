function [meanEn, confI] = getWindowsBandEnergy(band,varargin)
%  This function calculates the average energy of a signal sampled from
% windows of a specific band (all windows from a fixed distance from the cell edge)
%
% Note: Remove mean and linear trend before applying this function
%
% Usage: [meanEn, confI] = getWindowsBandEnergy()
%
% Input :
%       band  - biosensor signal band (# of windows X  # of time points) 
%
%       nBoot - # of boostrap samples to be used (default value 1000)  
%
%       alpha - alpha used to generate the bootstrap confidence intervals
%         (default value 0.05)
%
% Output:
%       meanEn - mean energy of the band
%
%       confI  - confidence interval    
%
% See also: removeMeanTrendNaN
%
%Marco Vilela, 2012

ip=inputParser;
ip.addRequired('band',@(x) iscell(x));
ip.addOptional('nBoot',1e3,@isscalar);
ip.addOptional('alpha',.05,@isscalar);

ip.parse(band,varargin{:});
nBoot = ip.Results.nBoot;
alpha = ip.Results.alpha;
 

avgEn = cell2mat( cellfun(@(x) mean(x.^2),band,'UniformOutput',0) );

opt = statset('UseParallel','never');
if matlabpool('size')
    opt = statset('UseParallel','always');
end

[confI,statEner] = bootci(nBoot,{@mean,avgEn},'alpha',alpha,...
          'type','bca','Options',opt);
      
meanEn = mean(statEner);