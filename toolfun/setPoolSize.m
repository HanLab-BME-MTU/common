function setPoolSize(poolSize)
%SETPOOLSIZE sets parallel pool (matlabpool) to specified size irrespective of current status
%
% setPoolSize
% setPoolSize(poolSize)
%
% Input: poolSize - the parallel pool (matlabpool) size (number of
% workers). If not specified / empty then the default size is used, unless
% a pool already exists in which case nothing is done.
% 
%
% Hunter Elliott
% 12/2015

if nargin < 1
    poolSize = [];
end

%Get current poolsize
poolobj = gcp('nocreate'); %Just check the current size
if isempty(poolobj)
    currSize = 0;
else
    currSize = poolobj.NumWorkers;
end

if isempty(poolSize)
    %This appears to be the only way to use the default poolsize setting,
    %and it throws an error if one already exists so we catch that.
    try
        parpool
    catch
    end
elseif currSize ~= poolSize
    
    %Delete the current pool (Doesn't seem to be possible to just
    %add/remove...)
    delete(poolobj)
        
    if poolSize > 0
        parpool(poolSize);
    end    
end
