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
    ip.addParameter('DisplayFunc',@defaultDisplayFunc,@(x) isa(x,'function_handle'));

    % Find first non-cell input to begin parsing parameters
    for inIdx = 1:length(varargin)
       if(~iscell(varargin{inIdx}))
            paramIdx = inIdx;
            break;
       elseif(inIdx == length(varargin))
           paramIdx = inIdx + 1;
       end
    end
    
    if(paramIdx == 1)
        error('parcellfun_progress:NotACell', ...
         'Input #2 expected to be a cell array.');
    end

    ip.parse(varargin{paramIdx:end});
    in = ip.Results;

    inSize = size(varargin{1});
    nRuns = numel(varargin{1});
    nCompleted = 0;
    out = cell(nRuns,nout);
    notComplete = true(nRuns,1);
    
    F = cellfun(@parfunc,varargin{1:paramIdx-1},'UniformOutput',false);
    F = [F{:}];
    cleanup = onCleanup(@()  cancel(F));
    
    startDT = datetime('now');
    nRunsTime = startDT;
    nProgressOut = in.DisplayFunc(nRuns,0,startDT,NaN,0);
        
    while(nCompleted < nRuns)
        completedIdx = [];
        currentOut = cell(1,nout);
        try
            [completedIdx,currentOut{:}] = fetchNext(F,in.UpdateInterval);
        catch err
%             disp(err);
            nCompleted = nCompleted + 1;
            completedIdx = find([F.Read] & notComplete',1);
            disp(completedIdx);
            [currentOut{:}] = in.ErrorHandler(err);
        end
        if(~isempty(completedIdx))
            notComplete(completedIdx) = false;
            out(completedIdx,:) = currentOut;
            nCompleted = nCompleted + 1;
            nRunsTime = datetime('now');
        end
        nProgressOut = in.DisplayFunc(nRuns,nCompleted,startDT,nRunsTime,nProgressOut);
    end
%     failed = ~cellfun(@isempty,{F.Error});
%     out(failed,:) = arrayfun(in.ErrorHandler,F(failed),'UniformOutput',false);
    
    % Convert 2D cell array into a 1D cell array of cells
    varargout = num2cell(out,1);
    if(in.UniformOutput)
        varargout = cellfun(@cell2mat,varargout,'UniformOutput',false);
        nArgOut = cellfun('prodofsize',varargout);
        assert(numel(nArgOut) == 1 && all(nArgOut(1) == nArgOut(2:end)), ...
            'parcellfun_progress:NonUniformOutput', ...
            'Non-uniform output with UniformOutput set to true');
    end
    % Make output reflect size of input
    varargout = cellfun(@(out) reshape(out,inSize),varargout,'UniformOutput',false);
        
    % Wrap func in parfeval for evaluation
    function F = parfunc(varargin)
        F = parfeval(func,nout,varargin{:});
    end

end
function nProgressOut = defaultDisplayFunc(nRuns,nCompleted,startDT,nRunsTime,nProgressOut)
bp = repmat('\b',1,nProgressOut);
    nowDT = datetime('now');
    sinceStart = nowDT - startDT;
    if(~nCompleted)
        estDuration = 'xx:xx:xx';
        remaining = 'xx:xx:xx';
    else
        estDuration = (nRunsTime-startDT)*nRuns/nCompleted;
        remaining = estDuration - sinceStart;
    end
    nProgressOut = fprintf([bp ...
        ' %g / %g\n' ...
        '  Estimated %s\n' ...
        '-   Elapsed %s\n' ...
        '---------------------\n' ...
        '= Remaining %s\n'],...
        nCompleted,nRuns, ...
        char(estDuration), ...
        char(sinceStart), ...
        char(remaining) ...
        ) - nProgressOut;
end

