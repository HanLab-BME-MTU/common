function bootstrapMoviesCorrelation(movieList,varargin)
% bootstrapMoviesCorrelation bootstrap the correlation over multiple movies
%
% SYNOPSIS bootstrapMoviesCorrelation(movieList,paramsIn)
%
% INPUT   
%   movieList - A MovieList object describing the list of movies to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson, Sep 2011
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieList', @(x) isa(x,'MovieList'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieList,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieList.getProcessIndex('CorrelationBootstrappingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieList.processes_)+1;
    movieList.addProcess(CorrelationBootstrappingProcess(movieList,...
        movieList.outputDirectory_));                                                                                                 
end

corrBootProc = movieList.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(corrBootProc,paramsIn);
p2.ProcessName = p.ProcessName;

%% --------------- Initialization ---------------%%

% Test the presence and valid correlationCalculationProcess
nMovies =numel(movieList.movies_);
corrProc=cell(nMovies,1);
for i=1:nMovies
    movieData = movieList.movies_{i};
    iCorrProc =movieData.getProcessIndex('CorrelationCalculationProcess',1,1);
    if isempty(iCorrProc)
        iCorrProc = numel(movieData.processes_)+1;
        movieData.addProcess(CorrelationCalculationProcess(movieData,...
            movieData.outputDirectory_));
    end
    corrProc{i} = movieData.processes_{iCorrProc};
    parseProcessParams(corrProc{i},p2);
    
    corrProc{i}.run();
end

%Parse input, store in parameter structure
p = parseProcessParams(corrAvgProc,paramsIn);


input = corrBootProc.getInput;
nInput=numel(input);
% Set up output files
inFilePaths=cell(nInput,nInput);
outFilePaths=cell(nInput,nInput);
for i=1:nInput
    for j=1:i-1
        inFilePaths{i,j} = corrProc.outFilePaths_{i,j};
        outFilePaths{i,j} = [p.OutputDirectory filesep 'correlation' ...
            input(i).name '_' input(j).name '.mat'];
    end
    inFilePaths{i,i} = corrProc.outFilePaths_{i,i};
    outFilePaths{i,i} = [p.OutputDirectory filesep 'autocorrelation' ...
        input(i).name '.mat'];
end
mkClrDir(p.OutputDirectory);
corrAvgProc.setInFilePaths(inFilePaths);
corrAvgProc.setOutFilePaths(outFilePaths);


%% --------------- Correlation averaging ---------------%%% 
disp('Starting averaging correlation...')

% Calculate cross-correlation for all data
for iInput1=1:nInput
    for iInput2=1:iInput1
        [corrFun,lags,bounds]=corrProc.loadChannelOutput(iInput1,iInput2,...
            'output',{'corrFun','lags','bounds'});
        
        %Filter out the windows for averaging
        p.windowsRange=6:20;
        corrFun = corrFun(:,p.windowsRange);
        lags = lags(:,p.windowsRange);
        bounds = bounds(:,p.windowsRange);
         
        % Perform averaging
        nSamples = sum(~arrayfun(@(i)all(isnan(corrFun(:,i))),1:size(corrFun,2)));
        s.avgCorrFun = nanmean(corrFun,2);
        s.meadianCorrFun = nanmedian(corrFun,2);
        s.stdCorrFun = nanstd(corrFun,0,2);
        s.steCorrFun = s.stdCorrFun/sqrt(nSamples);
        s.lags = nanmean(lags,2);
        s.avgBounds = nanmean(bounds,2);
        s.stdBounds = nanstd(bounds,0,2);
        save(outFilePaths{iInput1,iInput2},'-struct','s');
    end
end
disp('Finished averaging correlation...')
