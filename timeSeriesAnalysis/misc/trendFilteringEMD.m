function [trend,detrend] = trendFilteringEMD(TS,plotYes)
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
%       trend   - cell array with trends for each found scale
%       detrend - cell array (TS - trend)
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
ZCR       = zeroC(1:end-1)./zeroC(2:end);

ZCR(~isfinite(ZCR))           = [];
changeScale                   = find( (ZCR-limits(2) > 0 | ZCR -limits(1) < 0) );
changeScale(changeScale == 1) = [];

if isempty(changeScale)
    changeScale = size(imf,1)-1;
end

trend   = arrayfun(@(x) sum(imf(x+1:end,:),1),changeScale,'Unif',0);
detrend = cellfun(@(x) TS'-x,trend,'Unif',0);


if plotYes
    
    plotC = num2cell(hsv(numel(trend)),2);
    figure
    plot(TS,'-*')
    hold on
    cellfun(@(x,y) plot(x,'Color',y),trend,plotC,'Unif',0)
    
end

