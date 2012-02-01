function processMovieSignal(movieObject,varargin)
% PROCESSMOVIESIGNAL process time series from the sampled maps
%
% SYNOPSIS processMovieSignal(movieObject,paramsIn)
%
% INPUT   
%   movieObject - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson, Jan 2012

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObject', @(x) isa(x,'MovieObject'));
ip.addOptional('paramsIn',[], @isstruct);
ip.addParamValue('waitbar',[], @ishandle);
ip.parse(movieObject,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous  signal preprocessing process                                                                  
iProc = movieObject.getProcessIndex('SignalProcessingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieObject.processes_)+1;
    movieObject.addProcess(SignalProcessingProcess(movieObject,...
        movieObject.outputDirectory_));                                                                                                 
end

signalProc = movieObject.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(signalProc,paramsIn);

stack=dbstack;
if ~any(strcmp('Process.run',{stack(:).name}));
    signalProc.run();
    return;
end

%% --------------- Initialization ---------------%%
if ~isempty(ip.Results.waitbar)
    wtBar=ip.Results.waitbar;
    waitbar(0,ip.Results.waitbar,'Initializing...');
    wtBarArgs={'waitbar',wtBar};
elseif feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...');
    wtBarArgs={'waitbar',wtBar};
else
    wtBar=-1;
    wtBarArgs={};
end

% Set the waitbar title if applicable
if ishandle(wtBar),
    [~,movieName]=fileparts(movieObject.getPath);
    set(wtBar,'Name',movieName);
end

% Delegates signal processing to  individual movies if movie list 
if isa(movieObject,'MovieList')
    % Movie list function consists of two steps
    % 1- Delegation of process to individual movies (given by MovieIndex)
    % 2- Processing of concatenated data for all movies
    
    % Movie delegation
    movieParams=rmfield(p,{'MovieIndex','OutputDirectory'});
    nMovies = numel(p.MovieIndex);
    movieSignalProc=cell(nMovies,1);
    movieInput = cell(nMovies,1);
    
    
    for i =1:nMovies;
        iMovie = p.MovieIndex(i);
        fprintf(1,'Processing signal for movie %g/%g\n',i,nMovies);
        
        % Delegate  calculation for each movie of the list
        movieParams.SliceIndex=p.SliceIndex{iMovie};
        movieData = movieObject.movies_{iMovie};
        iProc = movieData.getProcessIndex('SignalProcessingProcess',1,0);
        if isempty(iProc)
            movieData.addProcess(SignalProcessingProcess(movieData,...
                movieData.outputDirectory_));
        end
        movieSignalProc{i} = movieData.processes_{end};
        parseProcessParams(movieSignalProc{i},movieParams);
        movieSignalProc{i}.run(wtBarArgs{:});
        movieInput{i}=movieSignalProc{i}.getInput;       
    end  
    
    % Check input is the same for all movies
    assert(all(cellfun(@(x) isequal(x,movieInput{1}),movieInput)));
    input=movieInput{1};
    nInput= numel(movieInput{1});
    
    % Check number of frames per movie and determine lag limits
    nFrames =  cellfun(@(x) x.nFrames_,movieObject.movies_);
    nFrames=max(nFrames);
    timeInterval=unique(cellfun(@(x) x.timeInterval_,movieObject.movies_));
    assert(numel(timeInterval)==1);
    
    % Load input for all  movies and concatenate them
    fprintf(1,'Processing signal for movie list.\n');
    [input.data] = deal({});
    [input.range] = deal({});
    disp('Using preprocessed signal from:');
    for i =1:nMovies;
        iMovie = p.MovieIndex(i);
        % Retrieve individual movie output
        movieInput{i} =loadInput(movieObject.movies_{iMovie},movieInput{i});
        for iInput=1:nInput
            localdata =  movieInput{i}(iInput).data;
            localrange =  movieInput{i}(iInput).range;
            for iBand=1:min(numel(input(iInput).data),numel(localdata))
                input(iInput).data{iBand} =vertcat(input(iInput).data{iBand},localdata{iBand});
                input(iInput).range{iBand} =vertcat(input(iInput).range{iBand},localrange{iBand});
            end
            for iBand=numel(input(iInput).data)+1:max(numel(input(iInput).data),numel(localdata))
                input(iInput).data{iBand} =localdata{iBand};
                input(iInput).range{iBand} =localrange{iBand};
            end
        end
    end
    inFilePaths = cellfun(@(x){x.paths},movieInput,'Unif',false);
    inFilePaths=vertcat(inFilePaths{:});
    signalProc.setInFilePaths(inFilePaths);
    p.SliceIndex = vertcat(p.SliceIndex{:});
else    
    
    % Initialization of MovieData object
    assert(isa(movieObject,'MovieData'));
    movieData=movieObject;
    
    % Retrieve time interval and number of frames
    timeInterval=movieData.timeInterval_;
    nFrames = movieData.nFrames_;
    assert(~isempty(timeInterval) && ~isempty(nFrames)) 
    
    % Set the waitbar title if applicable
    if ishandle(wtBar),
        [~,movieName]=fileparts(movieObject.getPath);
        set(wtBar,'Name',movieName);
    end
    
    input = signalProc.getInput;
    disp('Using preprocessed signal from:');
    input = loadInput(movieData,input);
    signalProc.setInFilePaths({input.paths});
end

%% Signal processing

% Set the waitbar title if applicable
if ishandle(wtBar),
    [~,movieName]=fileparts(movieObject.getPath);
    set(wtBar,'Name',movieName);
end

% Clear output directory and recreate fresh one
if exist(p.OutputDirectory,'dir'), rmdir(p.OutputDirectory,'s'); end
mkdir(p.OutputDirectory);

% Call individual processing tools
p.nFrames=nFrames;
p.timeInterval=timeInterval;
nInput= numel(input);
nTools = numel(p.processingTools);
outFilePaths =cell(nInput,nInput,nTools);
for i=1:nTools
    % Create output subdirectory for each tool
    allTools = signalProc.getProcessingTools;
    selectedTool = allTools(p.processingTools(i).type);
    toolParams=p.processingTools(i).parameters;
    toolParams.outputDir = fullfile(p.OutputDirectory,[num2str(i) '_' selectedTool.name]);
    mkClrDir(toolParams.outputDir);
    
    % Pass input, process parameters and tool specific parameters
    outFilePaths(:,:,i)=selectedTool.function(input,p,toolParams,wtBarArgs{:});
end
signalProc.setOutFilePaths(outFilePaths);

function input = loadInput(movieData,input)
% Load the output of the preprocessing for each input and add it to the
% input structure

nInput=numel(input);

% Test the presence and output validity of the signal preprocessing
iSignalPreproc =movieData.getProcessIndex('SignalPreprocessingProcess',1,1);
if isempty(iSignalPreproc)
    error([SignalPreprocessingProcess.getName ' has not yet been performed'...
        'on this movie! Please run first!!']);
end

% Check that there is a valid output
signalPreproc = movieData.processes_{iSignalPreproc};
preprocInput =signalPreproc.getInput;
preprocIndex=zeros(nInput,1);
for i=1:nInput
    index = find(arrayfun(@(x) isequal(input(i),x),preprocInput));
    assert(isscalar(index))
    preprocIndex(i) = index;
end
if ~signalPreproc.checkChannelOutput(preprocIndex)
    error(['Each time series must have been preprocessed !' ...
        'Please apply pre-processing to all time series before '...
        'running coherence calculation!'])
end

disp(signalPreproc.funParams_.OutputDirectory);

% Load input
for iInput=1:nInput
    input(iInput).paths = signalPreproc.outFilePaths_{1,preprocIndex(iInput)};
    [input(iInput).data,input(iInput).range] = signalPreproc.loadChannelOutput(preprocIndex(iInput));
end

