function [ varargout ] = parcellfun_progress( func, varargin )
%parcellfun_progress is a parallelized cellfun that also has
% a progress monitor built-in.
%
% See also cellfun, distributed.cellfun, parfor

% Mark Kittisopikul, December 2015
    nout = nargout;

    ip = inputParser;
    ip.addParameter('UniformOutput',true,@islogical);
    ip.addParameter('ErrorHandler',@rethrow,@(x) isa(x,'function_handle'));
    ip.addParameter('UpdateInterval',1,@(x) validateattributes(x,{'numeric'},{'scalar'}));

    % Find first non-cell input to begin parsing parameters
    for inIdx = 1:length(varargin)
       if(~iscell(varargin{inIdx}))
            break;
       end 
    end
    
    assert(inIdx ~= 1,'parcellfun_progress:NotACell','Input #2 expected to be a cell array.');

    ip.parse(varargin{inIdx:end});
    in = ip.Results;

    
    
    F = cellfun(@parfunc,varargin{:},'UniformOutput',false);
    F = [F{:}];
    
    nRuns = length(F);
    nCompleted = 0;
    out = cell(nRuns,nout);
    
    fprintf('Parallel Cell Function:');
    
    nProgressOut = fprintf(' %g / %g Estimating time ...\n',0,nRuns);
    backspaceChar = '\b';
    
    startDT = datetime('now');
    
    updateInterval = 1;
    
%     for idx = 1:nRuns
    while(nCompleted < nRuns)
        currentOut = cell(1,nout);
        [completedIdx,currentOut{:}] = fetchNext(F,updateInterval);
        if(~isempty(completedIdx))
            out(completedIdx,:) = currentOut;
            nCompleted = nCompleted + 1;
        end
        bp = repmat(backspaceChar,1,nProgressOut);
%         bp = '';
        nowDT = datetime('now');
        sinceStart = nowDT - startDT;
        estDuration = sinceStart*nRuns/nCompleted;
        remaining = estDuration - sinceStart;
        nProgressOut = fprintf([bp ' %g / %g\n Total %02g:%02g:%02g -  Elapsed %02g:%02g:%02g = Remaining %02g:%02g:%02g\n'],...
            nCompleted,nRuns, ...
            floor(hours(estDuration)), ...
            floor(minutes(estDuration)), ...
            ceil(seconds(estDuration)), ...
            floor(hours(sinceStart)), ...
            floor(minutes(sinceStart)), ...
            ceil(seconds(sinceStart)), ...
            floor(hours(remaining)), ...
            floor(minutes(remaining)), ...
            ceil(seconds(remaining)) ...
            ) - nProgressOut;
    end
    
    % Convert 2D cell array into a 1D cell array of cells
    varargout = num2cell(out,1);
    if(in.UniformOutput)
        varargout = cellfun(@cell2mat,varargout,'UniformOutput',false);
        nArgOut = cellfun('prodofsize',varargout);
        assert(numel(nArgOut) == 1 || all(nArgOut(1) == nArgOut(2:end)), ...
            'parcellfun_progress:NonUniformOutput','Non-uniform output with UniformOutput set to true');
    end
        

    function F = parfunc(varargin)
        F = parfeval(func,nout,varargin{:});
    end

end

