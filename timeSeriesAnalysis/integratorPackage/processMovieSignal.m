function processMovieSignal(movieObjectOrProcess,varargin)
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
%       ('FieldName' -> possible values)
%
%       ('MovieIndex' -> array of integers)  If movieObject is a movie
%       list, specifies the index of the movies to be preprocessed.

%
%       ('OutputDirectory' -> character string)  A character
%       string specifying the directory to save the processed output.
%       Each processing tool will be saved in a separate subdirectory.
%
%       ('ProcessName'-> Cell array of character strings) The name of the
%       sampling processes to retrieve the sampled output from. All output
%       of the process will be preprocessed.
%
%       ('tools'-> array of structures) The processing tools to
%       apply to the input time-series. Each tools should contain a type
%       field indexing the type of tools(see getTools static
%       method of SignalProcessingProcess) and a parameters field giving
%       the parameters to use when processing the time-series.

% Sebastien Besson, Jan 2012 (last modified Feb 2012

%% Input

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObjectOrProcess', @isProcessOrMovieObject);
ip.addOptional('paramsIn',[], @isstruct);
ip.addParamValue('waitbar',[], @ishandle);
ip.parse(movieObjectOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get MovieData object and Process
[movieObject, signalProc] = getOwnerAndProcess(movieObjectOrProcess,'SignalProcessingProcess',true);

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
    % Movie list processin consists of two steps
    % 1- Delegation of processing to individual movies (in MovieIndex)
    % 2- Processing of data concatenated for all movies
    
    % Movie delegation
    movieParams=rmfield(p,{'MovieIndex','OutputDirectory'});
    nMovies = numel(p.MovieIndex);
    movieSignalProc=cell(nMovies,1);
    movieInput = cell(nMovies,1);
    
    
    for i =1:nMovies;
        iMovie = p.MovieIndex(i);
        fprintf(1,'Processing signal for movie %g/%g\n',i,nMovies);
        
        % Delegate processing for each movie of the list
        movieData = movieObject.getMovies{iMovie};
        
        [~,movieSignalProc{i}] = getOwnerAndProcess(movieData,'SignalProcessingProcess',true);
        
        % Parse parameters and run the process
        parseProcessParams(movieSignalProc{i},movieParams);
        movieSignalProc{i}.run(wtBarArgs{:});
        movieInput{i}=movieSignalProc{i}.getInput;       
    end  
    
    % Check input is the same for all movies
    assert(all(cellfun(@(x) isequal(x,movieInput{1}),movieInput)));
    input=movieInput{1};
    nInput= numel(movieInput{1});
    
    % Check number of frames per movie and determine lag limits
    nFrames =  cellfun(@(x) x.nFrames_,movieObject.getMovies);
    nFrames=max(nFrames);
    timeInterval=unique(cellfun(@(x) x.timeInterval_,movieObject.getMovies));
    assert(numel(timeInterval)==1);
    
    % Load input for all  movies and concatenate them
    fprintf(1,'Processing signal for movie list.\n');
    [input.data] = deal({});
    [input.range] = deal({});
    disp('Using preprocessed signal from:');
    for i =1:nMovies;
        iMovie = p.MovieIndex(i);
        % Retrieve individual movie output
        movieInput{i} =loadInput(movieObject.getMovies{iMovie},movieInput{i});
        
        % Concatenate all data
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
    
    % Load preprocessing output
    input = signalProc.getInput;
    disp('Using preprocessed signal from:');
    input = loadInput(movieData,input);
    signalProc.setInFilePaths({input.paths});
end

%% Signal processing
% Set the waitbar title
if ishandle(wtBar),
    [~,movieName]=fileparts(movieObject.getPath);
    set(wtBar,'Name',movieName);
end

% Clear output directory to recreate fresh one
if exist(p.OutputDirectory,'dir'), rmdir(p.OutputDirectory,'s'); end
mkdir(p.OutputDirectory);

% Call individual processing tools

nInput= numel(input);
nTools = numel(p.tools);
outFilePaths =cell(nInput,nInput,nTools);
for i=1:nTools
    % Retrieve i-th tool function from the static method
    allTools = signalProc.getTools;
    selectedTool = allTools(p.tools(i).type);
    toolParams=p.tools(i).parameters;
    
    % Create output subdirectory for each tool
    toolParams.OutputDirectory = fullfile(p.OutputDirectory,[num2str(i) '_' selectedTool.name]);
    toolParams.nFrames=nFrames;
    toolParams.timeInterval=timeInterval;
    mkClrDir(toolParams.OutputDirectory);
    
    % Pass input, generic process parameters and tool specific parameters
    outFilePaths(:,:,i)=selectedTool.function(input,toolParams,wtBarArgs{:});
end
% Set the outputDirectory for the process
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
if ~signalPreproc.checkOutput(preprocIndex)
    error(['Each time series must have been preprocessed !' ...
        'Please apply pre-processing to all time series before '...
        'running coherence calculation!'])
end

disp(signalPreproc.funParams_.OutputDirectory);

% Load input
for iInput=1:nInput
    input(iInput).paths = signalPreproc.outFilePaths_{1,preprocIndex(iInput)};
    [input(iInput).data,input(iInput).range] = signalPreproc.loadOutput(preprocIndex(iInput),'output',{'data','range'});
end

