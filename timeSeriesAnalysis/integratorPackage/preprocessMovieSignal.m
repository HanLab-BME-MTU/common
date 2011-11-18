function preprocessMovieSignal(movieObject,varargin)
% preprocessMovieSignal preprocess time series from the protrusion and activity maps
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

% Marco Vilela, Sep 2011
% Sebastien Besson, Nov 2011
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObject', @(x) isa(x,'MovieObject'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieObject,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
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

% Delegates correlation processes to movies if object is a movieList 
if isa(movieObject,'MovieList')
    movieParams.ProcessName=p.ProcessName;
    for i =1:numel(p.MovieIndex);
        movieData = movieObject.movies_{i};
        iProc = movieData.getProcessIndex('SignalPreprocessingProcess',1,0);
        if isempty(iProc)
            iProc = numel(movieData.processes_)+1;
            movieData.addProcess(SignalPreprocessingProcess(movieData,...
                movieData.outputDirectory_));
        end
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

input = signalPreProc.getInput;
nInput=numel(input);
inFilePaths = cell(nInput,1);
inData = cell(nInput,1);
for i=1:nInput
    proc = movieObject.processes_{input(i).processIndex};
    if isempty(input(i).channelIndex)
        inFilePaths{i,:} = proc.outFilePaths_{1};
        inData{i} = proc.loadChannelOutput('output',input(i).var);
        inData{i} = reshape(inData{i},size(inData{i},1),1,size(inData{i},2));
        
    else
        inFilePaths{i,:} = proc.outFilePaths_{1,input(i).channelIndex};
        inData{i} = proc.loadChannelOutput(input(i).channelIndex,'output',input(i).var);
    end
end

% Set up output files
outFilePaths=cell(1,nInput);
for i=1:nInput
    outFilePaths{1,i} = [p.OutputDirectory filesep input(i).name '.mat'];
end
mkClrDir(p.OutputDirectory);
signalPreProc.setOutFilePaths(outFilePaths);

%% --------------- Correlation calculation ---------------%%% 

disp('Starting preprocessing signal...')

% Check input have the same size
allSizes =cellfun(@size,inData,'Unif',false);
nSlices = unique(cellfun(@(x) x(1),allSizes));
nBands = cellfun(@(x) x(2),allSizes);
nPoints = unique(cellfun(@(x) x(3),allSizes));
% nBands(2)=1;
assert(isscalar(nSlices));
assert(isscalar(nPoints));
%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;


logMsg = @(iInput,iBand) ['Please wait, preprocessing signal of ' input(iInput).name...
    ' for band ' num2str(iBand)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;

for iInput=1:nInput
    
    disp(logMsg(iInput,1));
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput,1)); end
    data= cell(nBands(iInput),1);
    range=cell(nBands(iInput),1);
    inputData = inData{iInput};
    for iBand=1:nBands(iInput)
        % Initialize data and range for the corresponding band
        data{iBand}=cell(nSlices,1);
        range{iBand}=cell(nSlices,1);
       
        % Get data and remove outliers
        rawData =squeeze(inputData(:,iBand,:));
        rawData(detectOutliers(rawData,p.kSigma)) = NaN;
        
        %% Percentage of NaN
        validSlices = (nPoints-sum(isnan(rawData),2))>=minP;
        
        [data{iBand}(validSlices) ,range{iBand}(validSlices)] = ...
            removeMeanTrendNaN(rawData(validSlices,:)');
        if ishandle(wtBar), waitbar(iBand/nBands(iInput),wtBar,logMsg(iInput,iBand)); end
    end
    save(outFilePaths{1,iInput},'data','range');  
end


disp('Finished preprocessing signal...')
if ishandle(wtBar), close(wtBar); end

