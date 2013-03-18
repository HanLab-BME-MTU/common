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

lags = -maxLag:maxLag;
nObs = length(x);
%xIn  = testSign(x);
%yIn  = testSign(y);

numSTD = 2;

bounds = [numSTD;-numSTD]/sqrt(nObs);


pVal = [];
switch corrT
    
    case 'Pearson'
        %NaN and outliers have no influence. Well, not too much.
        
               
        SX = flipud(buffer(x,nObs,nObs - 1));%delay x(t-n)
        SY = flipud(buffer(y,nObs,nObs - 1));%delay y(t-n)

        SX(maxLag+2:end,:) = [];
        SY(maxLag+2:end,:) = [];
        
         SX = SX + tril(nan(size(SX)));
         SY = SY + tril(nan(size(SY)));
         Y  = repmat(y,1,maxLag+1)+ tril(nan(size(SY)))';
         X  = repmat(x,1,maxLag+1)+ tril(nan(size(SX)))';
         
        lagX = num2cell(SX',1);
        lagY = num2cell(SY',1);
        Y    = num2cell(Y,1);
        X    = num2cell(X,1);
        
        posLag = cell2mat(cellfun(@(x,y) robustfit(x,y,'bisquare'),lagX,Y,'Unif',0));
        negLag = cell2mat(cellfun(@(x,y) robustfit(x,y,'bisquare'),X,lagY,'Unif',0));
        normalizationR = cell2mat(cellfun(@(x,y) mad(x,1)/mad(y,1),lagX,Y,'Unif',0));
        normalizationL = cell2mat(cellfun(@(x,y) mad(x,1)/mad(y,1),X,lagY,'Unif',0));
        
        CCL   = normalizationR.*posLag(2,:);
        CCR   = normalizationL.*negLag(2,:);
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

        