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
    ip.addParameter('ParallelPool',[]);

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

    % Parse input parameters
    ip.parse(varargin{paramIdx:end});
    in = ip.Results;
    
    % Initialize parallel pool
    if(isempty(in.ParallelPool))
        in.ParallelPool = gcp;
    end

    % Initialize variables
    % Size of cell input
    inSize = size(varargin{1});
    % Number of elements in cell input (number of parallel tasks)
    nRuns = numel(varargin{1});
    % Cell array to hold all output
    out = cell(nRuns,nout);
    % Logical array to mark which jobs have completed
    notComplete = true(nRuns,1);
    
    % Start parallel execution
    startDT = datetime('now');
    F = cellfun(@parfunc,varargin{1:paramIdx-1},'UniformOutput',false);
    F = [F{:}];
    
    % When this function completes, user quits, or exception occurs,
    % cleanup
    cleanup = onCleanup(@()  cancel(F));

    % Number of workers completed
    nCompleted = 0;
    % Time it took to complete nRuns
    nCompletedTime = startDT;
    
    % Number of bytes output by progress, for backspacing
    nProgressOut = in.DisplayFunc(nRuns,0,startDT,NaN,0);
    
    % We will not update the estimate until we have processed all finished
    % workers
    last_nCompleted = nCompleted;
    last_nRunsTime = nCompletedTime;
   
    
    % Main output checking loop
    while(nCompleted < nRuns)
        completedIdx = -1;
        % Process finished workers available at the moment
        while(~isempty(completedIdx))
            currentOut = cell(1,nout);
            try
                [completedIdx,currentOut{:}] = fetchNext(F,1);
            catch err
                completedIdx = find([F.Read] & notComplete',1);
                [currentOut{:}] = in.ErrorHandler(err);
            end
            if(~isempty(completedIdx))
                notComplete(completedIdx) = false;
                out(completedIdx,:) = currentOut;
                nCompleted = nCompleted + 1;
                nCompletedTime = datetime('now');
            end
            nProgressOut = in.DisplayFunc(nRuns,last_nCompleted,startDT,last_nRunsTime,nProgressOut);
        end
        % New estimate completed, update timers
        last_nCompleted = nCompleted;
        last_nRunsTime = nCompletedTime;
        nProgressOut = in.DisplayFunc(nRuns,nCompleted,startDT,nCompletedTime,nProgressOut);
    end
    
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
        F = parfeval(in.ParallelPool,func,nout,varargin{:});
    end

end
function nProgressOut = defaultDisplayFunc(nRuns,nCompleted,startDT,nCompletedTime,nProgressOut)
    bp = repmat('\b',1,nProgressOut);
    nowDT = datetime('now');
    sinceStart = nowDT - startDT;
    if(~nCompleted)
        estDuration = 'xx:xx:xx';
        remaining = 'xx:xx:xx';
    else
        estDuration = (nCompletedTime-startDT)*nRuns/nCompleted;
        remaining = estDuration - sinceStart;
    end
    if(nRuns ~= nCompleted)
        nProgressOut = fprintf([bp ...
            ' %g / %g = %3.0f%%\n' ...
            '  Estimated %s\n' ...
            '-   Elapsed %s\n' ...
            '---------------------\n' ...
            '= Remaining %s\n'],...
            nCompleted,nRuns,nCompleted/nRuns*100, ...
            char(estDuration), ...
            char(sinceStart), ...
            char(remaining) ...
            ) - nProgressOut;
    else
        nProgressOut = fprintf([bp ...
            ' %g / %g' ...
            ' Elapsed %s\n'], ...
            nCompleted,nRuns, ...
            char(sinceStart) ...
            ) - nProgressOut;
    end
end

