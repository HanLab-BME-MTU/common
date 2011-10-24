function [out,trend]=preWhitening(TS,method)
%This function removes any deterministic component from the input time
%series TS
%
%Synopsis:
%         [out,trend]=preWhitening(TS,method)  
%
%Input: TS     - matrix(# of observations,# of variables)
%       method - string - 'imf' or 'ar' (Default is 'imf')
%
%Output
%       out   - detrend time series
%       trend - well, this is the trend
%
%Reference
%Z. Wu, N. E. Huang, S. R. Long and C.-K. Peng, On the trend, detrending, and the
%variability of nonlinear and non-stationary time series, Proc. Natl. Acad. Sci. USA
%104 (2007) 14889?14894
%
% Marco Vilela, 2011

if nargin < 2
    method = 'imf';
end

[nobs,nvar] = size(TS);
max_order   = 8;
trend       = [];
out         = TS;
h           = inf(2,nvar);

for i=1:nvar
    
    h(1,i) = kpsstest(TS(:,i));%Ho is stationary
    h(2,i) = vratiotest(TS(:,i));%Ho is a random walk
    
    if or(h(1,i),~h(2,i))
        switch method
            case 'ar'
                ts       = iddata(TS(:,i),[],1);
                model    = arxstruc(ts,ts,[1:max_order]');
                bestM    = selstruc(model,'aic');
                finalM   = ar(ts,bestM(1));
                IR       = polydata(finalM);
                out(:,i) = filter(IR,1,TS(:,i));
            case 'imf'
                imf = empiricalModeDecomp( TS(:,i) ) ;
                
                for j=1:length(imf)/2
                    trend = sum( cat(1, imf{ end-j:end } ) )';
                    nTs   = TS(:,i) - trend;
                    h     = kpsstest(nTs);
                    if ~h
                        out(:,i) = nTs;
                        disp('Changed 1')
                        break;
                    end
                end
                
                if h
                    disp('Changed 2')
                    out(:,i) = preWhitening(TS(:,i),1);
                    h = kpsstest(out(:,i));
                end
        end
    end
end
out = out - repmat(mean(out),nobs,1);