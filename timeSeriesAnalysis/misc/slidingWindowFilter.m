function [varargout] = slidingWindowFilter(TS,winSize,operation)
%This function applies the filter's operation to TS sliding window of size
%winSize
%Usage:
%
%      [out] = slidingWindowFilter(TS,winSize,operation)
%
%Input:
%
%      TS        - vector or matrix contained the Time Series 
%      winSize   - length of the sliding window
%      operation - anonymous function with the filter kernel.  
%       
%       Ex: out = slidingWindowFilter(TS,10,@(x) mean(x,2))
%                 mean sliding filter
%
%           out = slidingWindowFilter(TS,10,@(x) var(x,[],2))
%                 mean sliding filter
%           
%       IMPORTANT - operations are done over the rows. That is the reason for applying the mean on the dimension 2
%Output:
%
%       out - filtered signal.Same size as TS
%
%Marco Vilela, 2012


if ~isa(operation,'function_handle')
    error('The filter kernel has to be a function handle');
end

[nObs,nVar]  = size(TS);

if winSize >= nObs 
    error('Window size is too large')
end

nOut   = nargout;
workTS = num2cell(TS,1);
bound  = floor(winSize/2) + 1;
kernel = cellfun(@(x) formatFilterInput(x,bound,winSize),workTS,'Unif',0);

%Filtering
if nVar == 1
    
    [varargout{1:nOut}] = operation(kernel{1});
    
else
    
    input               = cellfun(@(x) num2cell(x,2),kernel,'Unif',0);
    [varargout{1:nOut}] = cellfun(@(x,y) operation(x,y),input{1},input{2},'Unif',0);
    
end

end %End of main function

function H1 = formatFilterInput(TS,bound,winSize)

%Adding reflective boundary condition
TS = [flipud(TS(2:bound));TS;flipud(TS(end - bound + 1:end - 1))];

%Sliding Window of size winSize
H  = hankel(TS);
H1 = H(1:end-2*(bound-1),1:winSize);

end