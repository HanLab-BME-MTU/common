function setPoolSize(poolSize)
%SETPOOLSIZE sets parallel pool (matlabpool) to specified size irrespective of current status
% setPoolSize(poolSize)
%
% Input: poolSize - the parallel pool (matlabpool) size (number of workers)
% 
%
% Hunter Elliott
% 12/2015

%Get current poolsize
poolobj = gcp('nocreate'); %Just check the current size
if isempty(poolobj)
    currSize = 0;
else
    currSize = poolobj.NumWorkers;
end

if currSize ~= poolSize
    
    %Delete the current pool (Doesn't seem to be possible to just
    %add/remove...)
    delete(poolobj)
        
    if poolSize > 0
        parpool(poolSize);
    end    
end
