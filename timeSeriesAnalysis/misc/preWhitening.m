function [out,trend,imf] = preWhitening(TS,varargin)
%This function removes any deterministic component from the input time
%series TS
%
%Synopsis:
%         [out,trend,imf]=preWhitening(TS,varargin)  
%
%Input: TS          - matrix(# of observations,# of variables)
%       method      - string - 'imf' or 'ar' (Default is 'imf')
%       alpha       - alpha value for the stationary tests    
%       removeMean  - logical - true to remove the time series mean value
%
%Output
%       out   - detrend time series
%       trend - this is the sum of all derterministic components of the
%               signal
%       imf   - cell array with the intrinsic mode functions
%
%Reference
%Z. Wu, N. E. Huang, S. R. Long and C.-K. Peng, On the trend, detrending, and the
%variability of nonlinear and non-stationary time series, Proc. Natl. Acad. Sci. USA
%104 (2007) 14889-14894
%
%See also : removeMeanTrendNaN
%
% Marco Vilela, 2011

%%Parsing input
ip=inputParser;
ip.addRequired('TS',@isnumeric);
ip.addParamValue('method','imf',@ischar);
ip.addParamValue('alpha',0.01,@isnumeric);
ip.addParamValue('removeMean',false,@islogical);

ip.parse(TS,varargin{:})
method     = ip.Results.method;
alpha      = ip.Results.alpha;
removeMean = ip.Results.removeMean;

%% Initialization
[nVar,nObs] = size(TS);
maxOrder    = 8;%Upper bound on the arx model order - empirical
trend       = nan(nVar,nObs);
out         = TS;
imf         = [];
hypot       = inf(2,nVar);

%%
for iVar=1:nVar
    
    hypot(1,iVar) = kpsstest(TS(iVar,:));%Ho is stationary - 0 
    hypot(2,iVar) = vratiotest(TS(iVar,:),'alpha',alpha);%Ho is a random walk - 1
    
    if or(hypot(1,iVar),~hypot(2,iVar))
        switch method
            case 'ar'
                
                ts       = iddata(TS(iVar,:),[],1);
                model    = arxstruc(ts,ts,(1:maxOrder)');
                bestM    = selstruc(model,'aic');
                finalM   = ar(ts,bestM(1));
                IR       = polydata(finalM);
                
                out(iVar,:) = filter(IR,1,TS(iVar,:));
                
            case 'imf'
                
                imf  = emd(TS(iVar,:));
                %Testing each IMF
                nImf = size(imf,1);
                rW   = inf(1,nImf);
                sS   = inf(1,nImf);
                
                for j = 1:nImf
                    rW(j) = vratiotest( imf(j,:) );
                    sS(j) = kpsstest(   imf(j,:) );
                end
                %
                
                range =  find( ~rW | sS ) ;
                if ~isempty(range)
                    for j=numel(range):-1:1
                        trend = sum( imf( range(j:end), : ), 1 );
                        DTs   = TS(iVar,:) - trend;
                        hypot     = kpsstest(DTs);
                        if ~hypot
                            out(iVar,:) = DTs;
                            break;
                        end
                    end
                end
                
                if hypot
                    
                    out(iVar,:) = preWhitening(TS(iVar,:),'method','ar');
                    hypot           = kpsstest(out(iVar,:));
                    
                end
        end
    end
end
%%
if removeMean
    out = out - repmat(mean(out),1,nObs);
end