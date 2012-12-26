function [xCorr,bounds,lags,pVal] = nanCrossCorrelation(x,y,varargin)
%Inputs are column vectors
%Under Construction
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('x',@isvector);
ip.addRequired('y',@isvector);
ip.addParamValue('corrType','Pearson', @ischar)
ip.addParamValue('maxLag',0,@isscalar);
ip.addParamValue('local',numel(x)-1,@isscalar);

ip.parse(x,y,varargin{:});
local    = ip.Results.local;
maxLag   = ip.Results.maxLag;
corrT    = ip.Results.corrType;

x = x(:);
y = y(:);

lags = [-maxLag:maxLag];
nObs = length(x);
xIn  = testSign(x);
yIn  = testSign(y);
numSTD = 2;

bounds = [numSTD;-numSTD]/sqrt(nObs);


pVal = [];
switch corrT
    
    case 'Pearson'
        %NaN have no influence
        xIn(isnan(xIn)) = 0;
        yIn(isnan(yIn)) = 0;
        nObs            = min([sum(~isnan(xIn)) sum(~isnan(yIn))]);
        normalization   = ( (nObs - 1)*nanstd(x)*nanstd(y) )^-1;
        
        SX = flipud(buffer(xIn,nObs,nObs - 1));%delay x(t-n)
        SY = flipud(buffer(yIn,nObs,nObs - 1));%delay y(t-n)

        SX(maxLag+2:end,:) = [];
        SY(maxLag+2:end,:) = [];
        
        CCR   = normalization*xIn'*SY';
        CCL   = normalization*yIn'*SX';
        xCorr = [fliplr(CCR) CCL(2:end)];
        pVal  = pvalPearson('b',xCorr,nObs);
    case 'Kendall'
        
        [xCorr(:,1)] = modifiedKendallCorr(x,y,'local',local,'maxLag',maxLag);
        
    case 'Spearman'
        
end
xCorr = xCorr(:);

end %ENd of main function

function [TSout] = testSign(TS)
%Removes the mean value if the TS oscillates between positive and negative values
TSout = TS;
positiveTest = sum( TS(isfinite(TS)) ) == sum( TS(isfinite(TS))>0 );
negativeTest = sum( TS(isfinite(TS)) ) == sum( TS(isfinite(TS))<0 );
if ~xor(positiveTest,negativeTest)
    TSout = TS - nanmean(TS);
end

end %end testSign
    
function p = pvalPearson(tail, xCorr, nObs)
%
t = xCorr.*sqrt((nObs-2)./(1-xCorr.^2)); % 
switch tail
    case 'b'
        
        p = 2*tcdf(-abs(t),nObs-2);
        
    case 'r'
        
        p = tcdf(-t,nObs-2);
        
    case 'l'
        
        p = tcdf(t,nObs-2);
        
end

end

        