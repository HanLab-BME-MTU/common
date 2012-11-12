function [out] = slidingWindowFilter(TS,winSize,operation)
%This function applies the filter's operation to TS sliding window of size
%winSize
%Usage:
%      [out] = slidingWindowFilter(TS,winSize,operation)
%
%Input:
%      TS        - vector or matrix contained the Time Series 
%      winSize   - 
%      operation - anonymous function with the filter kernel.  
%
%Output:
%       out - filtered signal.Same size as TS
%Marco Vilela, 2012

out = [];

if ~isa(operation,'function_handle')
    error('The filter kernel has to be a function handle');
end

[nObs,nVar]  = size(TS);

if winSize >= nObs 
    error('Window size is too large')
end


workTS = num2cell(TS,1);
bound  = floor(winSize/2) + 1;
funOut = nargout(operation);
Hat    = cellfun(@(x) formatFilterInput(x,bound,winSize),workTS,'Unif',0);

%Filtering
if nVar == 1
    
    out   = operation(Hat{1});
    
else
    
    input = cellfun(@(x) num2cell(x,2),Hat,'Unif',0);
    out   = cell2mat(cellfun(@(x,y) funOut,input{1},input{2},'Unif',0));
    
end

end %End of main function

function H1 = formatFilterInput(TS,bound,winSize)

%Adding reflective boundary condition
TS = [flipud(TS(2:bound));TS;flipud(TS(end - bound + 1:end - 1))];
%Sliding Window of size winSize
H  = hankel(TS);
H1 = H(1:end-winSize,1:winSize);
end