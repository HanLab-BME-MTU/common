function outTS = getTimeSeriesTrend(TS,varargin)
% This function removes time series trend
%
%USAGE
%       outTS = getTimeSeriesTrend(TS,varargin)
%
%Input:
%       TS - Time Series (#ofVariables,#ofPoints)
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
plotYes  = ip.Results.plotYes;
trendT   = ip.Results.trendType;

[nVar,nObs] = size(TS);
trend       = TS;
dTS         = TS;

for iVar = 1:nVar
    
    if ismember(trendT,0)
        
        % Remove sample mean 
        trend(iVar,:)    = repmat( nanmean(TS(iVar,:)),nObs,1 );
        deltaFit(iVar,:) = norminv((1-alpha/2),0,1)*repmat( nanstd(TS(iVar,:)),nObs,1 )/sqrt(nObs);
        dTS(iVar,:)      = TS(iVar,:) - trend(iVar,:);
        
    elseif ismember(trendT,[1 2 3])
        
        switch trendT
            case 1
                fitFun = @(b,x)(b(1)*x + b(2));
                bInit  = [1 0]; %Initial guess for fit parameters.                
            case 2
                fitFun = @(b,x)(b(1)*exp(b(2)*x));
                bInit  = [1 0]; %Initial guess for fit parameters.
            case 3
                fitFun = @(b,x)(b(1)*exp(b(2)*x))+(b(3)*exp(b(4)*x));
                bInit  = [1 0 1 0]; %Initial guess for fit parameters.
        end
        
        fitOptions = statset('Robust','on','MaxIter',500,'Display','off');
        [bFit,resFit,~,covFit,mseFit] = nlinfit([1:nObs],TS(iVar,:),fitFun,bInit,fitOptions);
        %Get confidence intervals of fit and fit values
        [trend(iVar,:),deltaFit] = nlpredci(fitFun,1:nObs,bFit,resFit,'covar',covFit,'mse',mseFit);
        dTS(iVar,:)              = TS(iVar,:) - trend(iVar,:);

               
    elseif trendT == 5
        workTS = gapInterpolation(TS(iVar,:),1);
        % Remove all deterministic components
        [dTS(iVar,:),trend(iVar,:)] = preWhitening(workTS);
        
    end
    
    if plotYes
        
        figure
        plot(TS(iVar,:))
        hold on
        plot(trend(iVar,:),'r')
        plot(trend(iVar,:)+deltaFit,'r--')
        plot(trend(iVar,:)-deltaFit,'r--')
        
    end
    
    outTS(iVar).trend     = trend(iVar,:);
    outTS(iVar).detrendTS = dTS(iVar,:);
    
end
