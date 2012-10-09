function [out] = slidingWindowFilter(TS,winSize,operation)
%This function applies the filter's operation to TS sliding window of size
%winSize
%Usage:
%      [out] = slidingWindowFilter(TS,winSize,operation)
%
%Input:
%      TS        - vector Time Series 
%      winSize   - 
%      operation - anonymous function with the filter kernel.  
%
%Output:
%       out - filtered signal.Same size as TS
%Marco Vilel, 2012

out = [];

if winSize >= numel(TS) 
    error('Window size is too large')
end

TS    = TS(:);
%Adding reflective boundary condition
bound = floor(winSize/2) + 1;
TS    = [flipud(TS(2:bound));TS;flipud(TS(end - bound + 1:end - 1))];
%Filtering
H     = hankel(TS);
H1    = H(1:end-winSize + 1,1:winSize);
out   = operation(H1);