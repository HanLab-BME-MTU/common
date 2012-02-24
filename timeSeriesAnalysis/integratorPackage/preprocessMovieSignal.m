function preprocessMovieSignal(movieObject,varargin)
% PREPROCESSMOVIESIGNAL preprocess time series from the sampled maps
%
% SYNOPSIS preprocessMovieSignal(movieObject,paramsIn)
%
% INPUT   
%   movieObject - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('MovieIndex' -> array of integers)  If movieObject is a movie
%       list, specifies the index of the movies to be preprocessed.
%
%       ('OutputDirectory' -> character string) Optional. A character
%       string specifying the directory to save the preprocessed output.
%       Each time-series will be saved as an individual mat file.
%
%       ('ProcessName'-> Cell array of character strings) The name of the
%       sampling processes to retrieve the sampled output from. All output
%       of the process will be preprocessed.
%
%       ('kSigma'-> positive integer) The multiplier to use when removing
%       outliers from time-series. See detectOutliers.m for more info.
%
%       ('trendType'-> integer) A value specifying the type of detrending
%       to apply to the time-series. See removeMeanTrendNan for a list of
%       acceptable values.

% Marco Vilela, Sep 2011
% Sebastien Besson, Nov 2011 (last modified Feb 2012)

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObject', @(x) isa(x,'MovieObject'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieObject,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous  signal preprocessing process                                                                  
iProc = movieObject.getProcessIndex('SignalPreprocessingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieObject.processes_)+1;
    movieObject.addProcess(SignalPreprocessingProcess(movieObject,...
        movieObject.outputDirectory_));                                                                                                 
end

signalPreProc = movieObject.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(signalPreProc,paramsIn);

% If movie list, delegates signal preprocessing to individual movies
if isa(movieObject,'MovieList')
    movieParams=rmfield(p,{'MovieIndex','OutputDirectory'});
    nMovies =numel(p.MovieIndex);
    for i =1:nMovies;
        movieData = movieObject.movies_{i};
        fprintf(1,'Processing signal for movie %g/%g\n',i,nMovies);
        
        % Create movie process if empty
        iProc = movieData.getProcessIndex('SignalPreprocessingProcess',1,0);
        if isempty(iProc)
            iProc = numel(movieData.processes_)+1;
            movieData.addProcess(SignalPreprocessingProcess(movieData,...
                movieData.outputDirectory_));
        end
        
        % Parse parameters and run movie process
        signalPreProc = movieData.processes_{iProc};
        parseProcessParams(movieData.processes_{iProc},movieParams);
        signalPreProc.run();
    end    
    return;
end    

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    [~,movieName]=fileparts(movieObject.getPath);
    wtBar = waitbar(0,'Initializing...','Name',movieName);
else
    wtBar=-1;
end

% Load input files
input = signalPreProc.getInput;
nInput=numel(input);
inFilePaths = cell(nInput,1);
inData = cell(nInput,1);
disp('Using sampled output from:');
for i=1:nInput
    proc = movieObject.processes_{input(i).processIndex};
    if isempty(input(i).channelIndex)
        inFilePaths{i} = proc.outFilePaths_{1};
        inData{i} = proc.loadChannelOutput('output',input(i).var);
        inData{i} = reshape(inData{i},size(inData{i},1),1,size(inData{i},2));
    else
        % Load output for given channelIndex and outputIndex
        inFilePaths{i} = proc.outFilePaths_{1,input(i).channelIndex};
        inData{i} = proc.loadChannelOutput(input(i).channelIndex,input(i).outputIndex,'output',input(i).var);
    end
    disp(inFilePaths{i})
end
signalPreProc.setInFilePaths(inFilePaths);

% Set up output files
outFilePaths=cell(1,nInput);
for i=1:nInput
    outFilePaths{1,i} = [p.OutputDirectory filesep input(i).name '.mat'];
end
disp('Saving results under:');
disp(p.OutputDirectory);
mkClrDir(p.OutputDirectory);
signalPreProc.setOutFilePaths(outFilePaths);

%% --------------- Signal pre-processing ---------------%%% kSigma
disp('Starting preprocessing signal...')

% Check input have the same size
allSizes =cellfun(@size,inData,'Unif',false);
nSlices = unique(cellfun(@(x) x(1),allSizes));
nBands = cellfun(@(x) x(2),allSizes);
nPoints = unique(cellfun(@(x) x(3),allSizes));
assert(isscalar(nSlices) && isscalar(nPoints));

%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

% Create log messages
logMsg = @(iInput,iBand) ['Please wait, preprocessing signal of ' input(iInput).name];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nBandsTot = sum(nBands);

for iInput=1:nInput    
    % Show log
    disp(logMsg(iInput));
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
    
    % Initialize data and range ouptut
    data= cell(nBands(iInput),1);
    [data{:}]=deal(cell(nSlices,1));
    range = data;
    
    for iBand=1:nBands(iInput)
        % Get iBand data and remove outliers
        rawData =squeeze(inData{iInput}(:,iBand,:));
        rawData(detectOutliers(rawData,p.kSigma)) = NaN;
        
        % Check percentage of NaN
        validSlices = (nPoints-sum(isnan(rawData),2))>=minP;
        [data{iBand}(validSlices) ,range{iBand}(validSlices)] = ...
            removeMeanTrendNaN(rawData(validSlices,:)',p.trendType);   
                
        % Update waitbar
        if ishandle(wtBar), 
            tj=toc;
            nj = sum(nBands(1:iInput-1))+ iBand;
            waitbar(nj/nBandsTot,wtBar,sprintf([logMsg(iInput) timeMsg(tj*nBandsTot/nj-tj)]));
        end   
    end
    save(outFilePaths{1,iInput},'data','range');  
end

disp('Finished preprocessing signal...')
if ishandle(wtBar), close(wtBar); end