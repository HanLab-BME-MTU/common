function outTS = getTimeSeriesTrend(TS,varargin)
% This function removes time series trend
%
%USAGE
%       outTS = getTimeSeriesTrend(TS,varargin)
%
%Input:
%       TS - Time Series (#ofPoints,#ofVariables)
%
%       trendType:
%                  0 - remove only the mean value
%                  1 - remove linear trend
%                  2 - remove exponential trend
%                  3 - remove double exponential trend
%                  4 - remove nonlinear local trend
%                  5 - remove all determinitic component (trend + seasonal)
%
%Output:
%       outTS(iVariable).trend     - estimated trend
%       outTS(iVariable).detrendTS - detrended time series
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('TS',@(x) isnumeric(x));
ip.addOptional('alpha',.05,@isscalar);
ip.addOptional('nSurr',100,@isscalar);
ip.addOptional('minLength',30,@isscalar);
ip.addOptional('plotYes',0,@isscalar);
ip.addOptional('trendType',0,@isscalar);

ip.parse(TS,varargin{:});
alpha    = ip.Results.alpha;
nSurr    = ip.Results.nSurr;
minLen   = ip.Results.minLength;
plotYes  = ip.Results.plotYes;
trendT   = ip.Results.trendType;

[nObs,nVar] = size(TS);


for iVar = 1:nVar
    
    if ismember(trendType,[0 1])
        
        % Remove sample means or linear trend
        trend(:,iVar) = nanmean(TS(:,iVar));
        dTS(:,iVar)   = TS(:,iVar) - trend(:,iVar);
        
    elseif ismember(trendType,[1 2 3])
        
        switch trendType
            case 1
                fitFun = @(b,x)(b(1)*x + b(2));
                bInit = [1 0]; %Initial guess for fit parameters.                
            case 2
                fitFun = @(b,x)(b(1)*exp(b(2)*x));
                bInit = [1 0]; %Initial guess for fit parameters.
            case 3
                fitFun = @(b,x)(b(1)*exp(b(2)*x))+(b(3)*exp(b(4)*x));
                bInit = [1 0 1 0]; %Initial guess for fit parameters.
        end
        
        fitOptions = statset('Robust','on','MaxIter',500,'Display','off');
        [bFit,resFit,jacFit,covFit,mseFit] = nlinfit(1:nObs,TS(:,iVar),fitFun,bInit,fitOptions);
        %Get confidence intervals of fit and fit values
        [fitValues,deltaFit] = nlpredci(fitFun,1:nObs,bFit,resFit,'covar',covFit,'mse',mseFit);
        

        
    elseif trendType == 5
        
        % Remove all deterministic components
        [dTS(:,iVar),trend(:,iVar),imf(:,iVar)] = preWhitening(dTS{iVar}(interval{iVar}));
        
    end
    
    if plotYes
        
        figure
        plot(TS(:,iVar))
        hold on
        plot(trend(:,iVar),'r')
        
    end
    
    outTS(iVar).trend     = trend(:,iVar);
    outTS(iVar).detrendTS = dTS(:,iVar);
    
end
