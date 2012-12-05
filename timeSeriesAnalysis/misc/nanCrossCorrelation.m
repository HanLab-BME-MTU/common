function [xCorr,bounds,pVal] = nanCrossCorrelation(x,y,type,mLag,neighborhood)
%Inputs are column vectors
%Under Construction
%
%Marco Vilela, 2012

ip = inputParser;
ip.addRequired('x',@isvector);
ip.addRequired('y',@isvector);
ip.addParamValue('type','Pearson', @ischar)
ip.addParamValue('mLag',10,@isscalar);

x = x(:);
y = y(:);

nObs = length(x);
xIn  = testSign(x);
yIn  = testSign(y);
numSTD = 2;

bounds = [numSTD;-numSTD]/sqrt(nObs);


pVal = [];
switch type
    
    case 'Pearson'
        %NaN have no influence
        xIn(isnan(xIn)) = 0;
        yIn(isnan(yIn)) = 0;
        nObs            = min([sum(~isnan(xIn)) sum(~isnan(yIn))]);
        normalization   = ( (nObs - 1)*nanstd(x)*nanstd(y) )^-1;
        
        SX = flipud(buffer(xIn,nObs,nObs - 1));%delay x(t-n)
        SY = flipud(buffer(yIn,nObs,nObs - 1));%delay y(t-n)

        SX(mLag+2:end,:) = [];
        SY(mLag+2:end,:) = [];
        
        CCR   = normalization*xIn'*SY';
        CCL   = normalization*yIn'*SX';
        xCorr = [fliplr(CCR) CCL(2:end)];
        pVal  = pvalPearson('b',xCorr,nObs);
    case 'Kendall'
        
        for iLag = 1:mLag
            
            [CCR(iLag),rPval(iLag)]   = modifiedCorr(xIn(1:end-iLag+1),yIn(iLag:end),neighborhood,'type',type);
            [CCL(iLag),lPval(iLag)]   = modifiedCorr(xIn(iLag:end),yIn(1:end-iLag+1),neighborhood,'type',type);
            
        end
        xCorr = [fliplr(CCR) CCL(2:end)];
        pVal  = [fliplr(rPval) lPval(2:end)];
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

        