function alignMovieEvents(movieObject,varargin)
% alignMovieEvents align sampled maps using predetermine events
%
% SYNOPSIS alignMovieEvents(movieObject,paramsIn)
%
% INPUT   
%   movieObject - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Kwonmoo Lee, 2011
% Sebastien Besson, Dec 2011
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObject', @(x) isa(x,'MovieObject'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieObject,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieObject.getProcessIndex('EventAlignmentProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieObject.processes_)+1;
    movieObject.addProcess(EventAlignmentProcess(movieObject,...
        movieObject.outputDirectory_));                                                                                                 
end

eventProc = movieObject.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(eventProc,paramsIn);

% Delegates correlation processes to movies if object is a movieList 
if isa(movieObject,'MovieList')
    movieParams.ProcessName=p.ProcessName;
    for i =1:numel(p.MovieIndex);
        movieData = movieObject.movies_{p.MovieIndex(i)};
        iProc = movieData.getProcessIndex('EventAlignmentProcess',1,0);
        if isempty(iProc)
            iProc = numel(movieData.processes_)+1;
            movieData.addProcess(EventAlignmentProcess(movieData,...
                movieData.outputDirectory_));
        end
        movieEventProc = movieData.processes_{iProc};
        parseProcessParams(movieData.processes_{iProc},movieParams);
        movieEventProc.run();
    end   
    
    % Here comes the bootstrapping stuff
    return;
end    



%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    [~,movieName]=fileparts(movieObject.getPath);
    wtBar = waitbar(0,'Initializing...','Name',movieName);
else
    wtBar=-1;
end

% 
% % Test the presence and output validity of the speckle detection process
% iSignalPreproc =movieObject.getProcessIndex('SignalPreprocessingProcess',1,1);     
% if isempty(iSignalPreproc)
%     error([SignalPreprocessingProcess.getName ' has not yet been performed'...
%     'on this movie! Please run first!!']);
% end        
% 
% % %Check that there is a valid output
% signalPreproc = movieObject.processes_{iSignalPreproc};
% preprocInput =signalPreproc.getInput;
% preprocIndex=zeros(nInput,1);
% for i=1:nInput
%     index = find(arrayfun(@(x) isequal(input(i),x),preprocInput));
%     assert(isscalar(index))
%     preprocIndex(i) = index;
% end
% if ~signalPreproc.checkChannelOutput(preprocIndex)
%     error(['Each time series must have been preprocessesd !' ...
%         'Please apply pre-processing to all time series before '...
%         'running correlation calculatino!'])
% end

% Load input files
input = eventProc.getInput;
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
eventProc.setInFilePaths(inFilePaths);


% Set up output files
outFilePaths=cell(nInput,nInput);
for i=1:nInput
    outFilePaths{1,i} = [p.OutputDirectory filesep 'aligned_' ...
            input(i).name '.mat'];
end
mkClrDir(p.OutputDirectory);
eventProc.setOutFilePaths(outFilePaths);

%% --------------- Correlation calculation ---------------%%% 
disp('Starting aligning events...')

%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
nFrames = movieObject.nFrames_;
nLagsMax =round(nFrames/4);

logMsg = @(i) ['Please wait, aligning ' input(i).name];
tic;

% Retrieve event times master process events
nWindows=size(inData{1},1);
eventTimes = NaN(nWindows+1,1);

eventData =squeeze(inData{1}(:,1,:));
[~,time]=max(eventData(:,nLagsMax:nFrames-nLagsMax),[],2);
validTimes= (time>1) & (time<(nFrames-2*nLagsMax+1));
eventTimes(validTimes)=time(validTimes)+nLagsMax-1;



% Align maps with regards to these events
alignedData = cell(nInput,1);
for iInput=1:nInput
    inputSize = size(inData{iInput});
    alignedData{iInput} =  nan(inputSize(1),inputSize(2),2*nFrames+1);
end 
if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end

disp(logMsg(iInput));
for iWindow=find(~isnan(eventTimes))'
    range = (nFrames+2-eventTimes(iWindow)):(2*nFrames+1-eventTimes(iWindow));
    for iInput=1:nInput
        alignedData{iInput}(iWindow,:,range)=inData{iInput}(iWindow,:,:);
    end
%     save(outFilePaths,'alignedData');  
end
    save(outFilePaths,'alignedData');  
disp('Finished aligning events...')
if ishandle(wtBar), close(wtBar); end

