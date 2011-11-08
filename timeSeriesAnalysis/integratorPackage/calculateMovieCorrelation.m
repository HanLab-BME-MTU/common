function calculateMovieCorrelation(movieObject,varargin)
% calculateMovieCorrelation calculate the autocorrelation and cross-correlation
% between the protrusion and activity maps
%
% SYNOPSIS calculateMovieCorrelation(movieObject,paramsIn)
%
% INPUT   
%   movieObject - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Marco Vilela, Sep 2011
% Sebastien Besson, Sep 2011
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObject', @(x) isa(x,'MovieObject'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieObject,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieObject.getProcessIndex('CorrelationCalculationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieObject.processes_)+1;
    movieObject.addProcess(CorrelationCalculationProcess(movieObject,...
        movieObject.outputDirectory_));                                                                                                 
end

corrProc = movieObject.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(corrProc,paramsIn);

% Delegates correlation processes to movies if object is a movieList 
if isa(movieObject,'MovieList')
    movieParams.ProcessName=p.ProcessName;
    for i =1:numel(p.MovieIndex);
        movieData = movieObject.movies_{i};
        iProc = movieData.getProcessIndex('CorrelationCalculationProcess',1,0);
        if isempty(iProc)
            iProc = numel(movieData.processes_)+1;
            movieData.addProcess(CorrelationCalculationProcess(movieData,...
                movieData.outputDirectory_));
        end
        corrProc = movieData.processes_{iProc};
        parseProcessParams(movieData.processes_{iProc},movieParams);
        corrProc.run();
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

input = corrProc.getInput;
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
outFilePaths=cell(nInput,nInput);
for i=1:nInput
    for j=1:i-1
        outFilePaths{i,j} = [p.OutputDirectory filesep 'correlation' ...
            input(i).name '_' input(j).name '.mat'];
    end
    outFilePaths{i,i} = [p.OutputDirectory filesep 'autocorrelation' ...
        input(i).name '.mat'];
end
mkClrDir(p.OutputDirectory);
corrProc.setOutFilePaths(outFilePaths);

%% --------------- Correlation calculation ---------------%%% 

disp('Starting calculating correlation...')

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
nLagsMax =round(nPoints/4);
data=cell(nInput,1);
range=cell(nInput,1);

logMsg = @(i) ['Please wait, calculating ' input(i).name ' autocorrelation'];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;

for iInput=1:nInput
    disp(logMsg(iInput));
    % Initialize autocorrelation function and bounds
    corrFun = nan(nLagsMax+1,nSlices,nBands(iInput));
    lags = nan(nLagsMax+1,nSlices,nBands(iInput));
    bounds = nan(2,nSlices,nBands(iInput));
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
    
    for iBand=1:nBands(iInput)
        % Initialize data and range for the corresponding band
        data{iInput}{iBand}=cell(nSlices,1);
        range{iInput}{iBand}=cell(nSlices,1);
       
        % Get data and remove outliers
        rawData =squeeze(inData{iInput}(:,iBand,:));
        rawData(detectOutliers(rawData,5)) = NaN;
        
        %% Percentage of NaN
        validSlices = (nPoints-sum(isnan(rawData),2))>=minP;
        
        [data{iInput}{iBand}(validSlices) ,range{iInput}{iBand}(validSlices)] = ...
            removeMeanTrendNaN(rawData(validSlices,:)');
        
        % Calculate autocorelation
        %         for iSlice=1:nSlices
        %             if length(data{iInput}{iBand}{iSlice}) >= minP
        
        for iSlice=find(validSlices)'
            nLags = round(length(data{iInput}{iBand}{iSlice})/4);
            [corrFun(1:nLags+1,iSlice,iBand),lags(1:nLags+1,iSlice,iBand),...
                bounds(:,iSlice,iBand)] = autocorr(data{iInput}{iBand}{iSlice},nLags);
%             end
        end
        if ishandle(wtBar), waitbar(iBand/nBands(iInput),wtBar,logMsg(iInput)); end
    end
    lags =lags*movieObject.timeInterval_; %#ok<NASGU>
    save(outFilePaths{iInput,iInput},'corrFun','bounds','lags');  
end

logMsg = @(i,j) ['Please wait, calculating ' input(i).name '/'...
    input(j).name ' cross-correlation'];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;

% Calculate cross-correlation for all data
for iInput1=1:nInput
    for iInput2=1:iInput1-1
    % Initialize cross-correlation function and bounds
        corrFun = nan(2*nLagsMax+1,nSlices,nBands(iInput1),nBands(iInput2));
        bounds  = nan(2,nSlices,nBands(iInput1),nBands(iInput2));
        lags  = nan(2*nLagsMax+1,nSlices,nBands(iInput1),nBands(iInput2));
        
        disp(logMsg(iInput1,iInput2));
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput1,iInput2)); end
        for iBand1=1:nBands(iInput1)
            for iBand2=1:nBands(iInput2)
                for iSlice=1:nSlices
                    [~,range1,range2] = intersect(range{iInput1}{iBand1}{iSlice},range{iInput2}{iBand2}{iSlice});
                    ccL               = length(range1);
                    if ccL >= minP
                        nLags = round(ccL/4);
                        [corrFun(1:2*nLags+1,iSlice,iBand1,iBand2),lags(1:2*nLags+1,iSlice,iBand1,iBand2),...
                            bounds(:,iSlice,iBand1,iBand2)] =...
                            crosscorr(data{iInput1}{iBand1}{iSlice}(range1),data{iInput2}{iBand2}{iSlice}(range2),nLags);
                    end
                end
            end
        end
        lags=lags*movieObject.timeInterval_; %#ok<NASGU>
        save(outFilePaths{iInput1,iInput2},'corrFun','bounds','lags');
    end
end
disp('Finished calculating correlation...')
if ishandle(wtBar), close(wtBar); end

