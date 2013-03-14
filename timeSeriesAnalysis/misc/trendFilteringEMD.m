function out = trendFilteringEMD(TS,plotYes)
%This function is an implementation of the trend filtering algorithm proposed in the reference
%
%Usage:  
%         out = trendFilteringEMD(TS,plotYes)
%
%Input:
%       TS       - time series vector
%       plotYes  - logical 
%
%Output
%       out - structure with fields: .detrend and .trend
%
%See Also: getTimeSeriesTrend
%
%Ref: Trend Filtering via Empirical Mode Decompositions. Moghtaderi at al. Computational Statistics and Data Analysis.
%Vol 58, Pag 114-126
%
%Marco Vilela
%2013

if nargin < 2
    plotYes = false;
end


%Empirical limits derived from zero-cross ratio IMF's obtained from 1000 simulations of fractional brownian motion with 1000 time points. See ref.
limits    = [1.93 2.4];

TS        = TS(:);
imf       = emd(TS);
zeroC     = getIMFzeroCrossing(imf);
pureTrend = imf(zeroC == 0,:);
ZCR       = zeroC(1:end-1)./zeroC(2:end);

ZCR(~isfinite(ZCR)) = [];

imfN = find(ZCR-limits(2) > 0);

slowScale = max(imfN);

if slowScale <= size(imf,1)
    trend     = sum(imf(slowScale+1:end,:),1)';
elseif ~isempty(pureTrend)
    trend     = pureTrend';
else
    trend     = nan(size(TS));
end

if plotYes
    
    figure
    plot(TS,'*')
    hold on
    plot(trend,'r')
    
end


out.dTS     = (TS - trend)';
out.trend   = trend';

