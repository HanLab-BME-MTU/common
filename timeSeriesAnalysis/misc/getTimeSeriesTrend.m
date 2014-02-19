function [outTS,exclude] = getTimeSeriesTrend(TS,varargin)
% This function removes time series trend.
% Accepts NaN (Of course, more NaN the time series has crappier the result is)
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
%                  4 - remove nonlinear local trend (trendFilteringEMD) - does not work for time series with NaN
%                  5 - remove all determinitic component and spits out a stationary signal
%                      the trend in this case is a smoothed version of the input signal - does not work for time series with NaN
%
%Output:
%       outTS(iVariable).trend     - estimated trend
%       outTS(iVariable).detrendTS - detrended time series
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('TS',@(x) isnumeric(x));
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('nSurr',100,@isscalar);
ip.addParamValue('plotYes',0,@isscalar);
ip.addParamValue('trendType',0,@isscalar);
ip.addParamValue('minLength',5,@isscalar);%minimum length. Ill-posed otherwise
ip.parse(TS,varargin{:});
alpha    = ip.Results.alpha;
plotYes  = ip.Results.plotYes;
trendT   = ip.Results.trendType;
minLen   = ip.Results.minLength;

%Initialization
[nVar,nObs] = size(TS);
trend       = TS;
dTS         = TS;
deltaFit    = nan(nVar,nObs);
iEx         = 1;
exclude     = nan(1,nObs);

for iVar = 1:nVar
    
    outTS.dTS(iVar,:)   = TS(iVar,:);
    outTS.trend(iVar,:) = nan(1,nObs);
    
    if sum(isfinite(TS(iVar,:))) > minLen
        
        if ismember(trendT,0)
            
            % Remove sample mean
            trend(iVar,:)    = repmat( nanmean(TS(iVar,:)),nObs,1 );
            deltaFit(iVar,:) = norminv((1-alpha/2),0,1)*repmat( nanstd(TS(iVar,:)),nObs,1 )/sqrt(nObs);
            dTS(iVar,:)      = TS(iVar,:) - trend(iVar,:);
            
        elseif ismember(trendT,[1 2 3])
            
            switch trendT
                case 1
                    fitFun = @(b,x)(b(1)*x + b(2));
                    bInit  = [nanmean(outTS.dTS(iVar,:)) randn(1)]; %Initial guess for fit parameters.
                case 2
                    fitFun = @(b,x)(b(1)*exp(b(2)*x));
                    bInit  = [nanmean(outTS.dTS(iVar,:)) randn(1)]; %Initial guess for fit parameters.
                case 3
                    fitFun = @(b,x)(b(1)*exp(b(2)*x))+(b(3)*exp(b(4)*x));
                    bInit  = rand(1,4); %Initial guess for fit parameters.
            end
            
            fitOptions = statset('Robust','on','RobustWgtFun','welsch','MaxIter',500,'Display','off');
            [bFit,resFit,~,covFit,mseFit] = nlinfit(1:nObs,TS(iVar,:),fitFun,bInit,fitOptions);
            
            %Get confidence intervals of fit and fit values
            [outTS.trend(iVar,:),deltaFit] = nlpredci(fitFun,1:nObs,bFit,resFit,'covar',covFit,'mse',mseFit);
            outTS.dTS(iVar,:)              = TS(iVar,:) - outTS.trend(iVar,:);
            
        elseif ismember(trendT,[4 5])
            
            
            nanTime      = isnan(TS(iVar,:));
            workTS       = gapInterpolation(TS(iVar,:),nObs);
            nanTimePost  = isnan(workTS);
            
            if trendT == 4
                
                [aux1,aux2]  = trendFilteringEMD(workTS(~nanTimePost));
                deTrend      = aux2{end};
                trendC       = aux1{end};
                
            else
                
                [deTrend,trendC] = preWhitening(workTS(~nanTimePost));
                
            end
            
            if isempty(find(nanTime,1))
                
                outTS.dTS(iVar,:)   = deTrend;
                outTS.trend(iVar,:) = trendC;
                
            else
                
                    
                newTime = find(~nanTime);    
                if nanTimePost(1) == 1
                    
                    [~,bLength] = findBlock(find(nanTimePost),1);
                    newTime     = find(~nanTime) - bLength(1);
                    
                end
                outTS.dTS(iVar,~nanTime)   = deTrend(newTime);
                outTS.trend(iVar,~nanTime) = trendC(newTime);
                
            end
            
        end
        
        if plotYes
            
            figure
            plot(TS(iVar,:))
            hold on
            plot(outTS.trend(iVar,:),'r')
            plot(outTS.trend(iVar,:)+deltaFit,'r--')
            plot(outTS.trend(iVar,:)-deltaFit,'r--')
            
        end
    else
        exclude(iEx) = iVar;
        iEx          = iEx + 1;
    end
    
end
exclude(isnan(exclude)) = [];