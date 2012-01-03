function calculateMovieCorrelation(movieObject,varargin)
% calculateMovieCorrelation calculate the correlation of temporal maps
%
% SYNOPSIS calculateMovieCorrelation(movieObject,paramsIn)
%
% INPUT   
%   movieObject - A movie object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       correlation function to.
%       If not input, the masks will be saved to the same directory as the
%       movieObject, in a sub-directory called "correlation"
%
%       ('MovieIndex' -> Positive integer scalar or vector)
%       Optional. For movie list, the integer index of the movie(s) to 
%       correlate. If not input, all movies in the list will be correlated. 
%
%       ('ProcessName' -> cell array of strings) Optional.
%       Names of processes to be correlated. The names should correspond to 
%       existing process classes.
% 
%       ('BandMin' -> Positive scalar)
%       The value of the minimal band of windows to be correlated.
%       Optional. Default is 0 (no jump suppression)
%
%       ('BandMin' -> Positive scalar)
%       The value of the maximal band of windows to be correlated.
%       Optional. Default is set to number of bands if input is a movie or 
%       the minimal value of the number of bands for all movies if the
%       input is a movie list
%
%       ('Slice Index' -> Logical array or cell array of logical arrays) 
%       A logical array of length equal to the number of window slices.
%       True indicates the slice should be included in the correlation. If
%       the input is a movie list, this should be a cell array of logical
%       array where the ith element of the cell array corresponds to the
%       slice index of the ith movie.
% 
%       ('nBoot' -> Positive integer)
%       The number of bootstrap data  samples to use during correlation
%       bootstrapping. Default is 1000.
%
%       ('alpha' -> Positive scalar)
%       The alpha used to generate the bootstrap confidence intervals
%       Default value is 0.05.
%
%

% Marco Vilela, Sep 2011
% Sebastien Besson, Sep 2011 (last modified Jan 2012)
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObject', @(x) isa(x,'MovieObject'));
ip.addOptional('paramsIn',[], @isstruct);
ip.addParamValue('waitbar',[], @ishandle);
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

stack=dbstack;
if ~any(strcmp('Process.run',{stack(:).name}));
    corrProc.run();
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

% Delegates correlation processes to movies if object is a movieList 
if isa(movieObject,'MovieList')
    movieParams=rmfield(p,{'MovieIndex','OutputDirectory'});
    nMovies= numel(p.MovieIndex);
    for i =1:nMovies;
        fprintf(1,'Calculating correlation for movie %g/%g\n',i,nMovies);
        
        % Delegate correlation calculation for each movie of the list
        movieParams.SliceIndex=p.SliceIndex{p.MovieIndex(i)};
        movieData = movieObject.movies_{p.MovieIndex(i)};
        iProc = movieData.getProcessIndex('CorrelationCalculationProcess',1,0);
        if isempty(iProc)
            iProc = numel(movieData.processes_)+1;
            movieData.addProcess(CorrelationCalculationProcess(movieData,...
                movieData.outputDirectory_));
        end
        corrProc = movieData.processes_{iProc};
        parseProcessParams(movieData.processes_{iProc},movieParams);
        corrProc.run(wtBarArgs{:});
    end  
    
    % Calls the movie list correlation bootstrapping method
    bootstrapMovieListCorrelation(movieObject,wtBarArgs{:});
    if ishandle(wtBar)&&isempty(ip.Results.waitbar), close(wtBar); end
    return;
end    

% This part should be only executed for MovieData objects
assert(isa(movieObject,'MovieData'));
movieData=movieObject;
assert(~isempty(movieObject.timeInterval_));

% Set the waitbar title if applicable
if ishandle(wtBar), 
    [~,movieName]=fileparts(movieObject.getPath);
    set(wtBar,'Name',movieName);
end

input = corrProc.getInput;
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
    error(['Each time series must have been preprocessesd !' ...
        'Please apply pre-processing to all time series before '...
        'running correlation calculatino!'])
end

% Load input
inFilePaths = cell(nInput,1);
data = cell(nInput,1);
range = cell(nInput,1);
for iInput=1:nInput
    inFilePaths{iInput,1} = signalPreproc.outFilePaths_{1,preprocIndex(iInput)};
    [data{iInput},range{iInput}] = signalPreproc.loadChannelOutput(preprocIndex(iInput));
end
disp('Using preprocessed signal from:');
disp(signalPreproc.funParams_.OutputDirectory);

corrProc.setInFilePaths(inFilePaths);

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
disp('Results will be saved under:')
disp(p.OutputDirectory);
mkClrDir(p.OutputDirectory);
corrProc.setOutFilePaths(outFilePaths);

%% --------------- Correlation calculation ---------------%%% 

%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nLagsMax =round(movieData.nFrames_/4);
nBands =cellfun(@numel,data);
nSlices = numel(data{1}{1});

logMsg = @(i) ['Please wait, calculating ' input(i).name ' autocorrelation'];

% Calculate autocorrelation
lags =(0:nLagsMax)'*movieData.timeInterval_; %#ok<NASGU>
for iInput=1:nInput
    disp(logMsg(iInput));
    
    % Initialize autocorrelation function and bounds
    corrFun = NaN(nLagsMax+1,nSlices,nBands(iInput));
    bounds = NaN(2,nSlices,nBands(iInput));
    bootstrapCorrFun=NaN(nLagsMax+1,nBands(iInput));
    bootstrapBounds=NaN(2,nLagsMax+1,nBands(iInput));
        
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
    
    for iBand=p.BandMin:min(nBands(iInput),p.BandMax);
        
        % Get number of timepoints and prune out slices
        nTimepoints=cellfun(@length,data{iInput}{iBand});
        validSlices =nTimepoints >=minP & p.SliceIndex;
        
        % Calculate raw auto-correlation
        for iSlice=find(validSlices)'
            nLags = round(length(data{iInput}{iBand}{iSlice})/4);
            [corrFun(1:nLags+1,iSlice,iBand),~,bounds(:,iSlice,iBand)] = ...
                autocorr(data{iInput}{iBand}{iSlice},nLags);
        end
        
        % Bootstrap valid autocorrelation functions
        validCorrFunSlices = sum(isnan(corrFun(:,:,iBand)),1)==0;
        if sum(validCorrFunSlices)>2
            [meanCC,CI] = correlationBootstrap(corrFun(:,validCorrFunSlices,iBand),...
                bounds(1,validCorrFunSlices,iBand),p.nBoot,p.alpha);
            bootstrapCorrFun(:,iBand)=meanCC;
            bootstrapBounds(:,:,iBand)=CI;
        end
        
        if ishandle(wtBar), waitbar(iBand/nBands(iInput),wtBar); end
    end
    
    save(outFilePaths{iInput,iInput},'corrFun','bounds','lags',...
        'bootstrapCorrFun','bootstrapBounds');  
end

logMsg = @(i,j) ['Please wait, calculating ' input(i).name '/'...
    input(j).name ' cross-correlation'];

% Calculate cross-correlation
lags =(-nLagsMax:nLagsMax)'*movieData.timeInterval_; %#ok<NASGU>
for iInput1=1:nInput
    for iInput2=1:iInput1-1
        disp(logMsg(iInput1,iInput2));
        
        % Initialize cross-correlation function and bounds
        corrFun = NaN(2*nLagsMax+1,nSlices,nBands(iInput1),nBands(iInput2));
        bounds  = NaN(2,nSlices,nBands(iInput1),nBands(iInput2));
        bootstrapCorrFun=NaN(2*nLagsMax+1,nBands(iInput1),nBands(iInput2));
        bootstrapBounds=NaN(2,2*nLagsMax+1,nBands(iInput1),nBands(iInput2));
        
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput1,iInput2)); end
        
        % Loop over bands and window slices
        bands1=p.BandMin:min(nBands(iInput1),p.BandMax);
        for i1=1:numel(bands1)
            iBand1=bands1(i1);
            for iBand2=p.BandMin:min(nBands(iInput2),p.BandMax)
               
                % Find valid range and test minimum number of timepoints
                nTimepoints = cellfun(@(x,y) length(intersect(x,y)),range{iInput2}{iBand2},...
                    range{iInput1}{iBand1});
                validSlices = nTimepoints>=minP & p.SliceIndex;
                
                % Calculate raw cross-correlation
                for iSlice=find(validSlices)'
                    % Retrieve number of lags from range intersection
                    [~,range1,range2] = intersect(range{iInput1}{iBand1}{iSlice},range{iInput2}{iBand2}{iSlice});
                    nLags = round(length(range1)/4);
                    [corrFun(1:2*nLags+1,iSlice,iBand1,iBand2),~,bounds(:,iSlice,iBand1,iBand2)] =...
                        crosscorr(data{iInput1}{iBand1}{iSlice}(range1),data{iInput2}{iBand2}{iSlice}(range2),nLags);
                end
                
                % Bootstrap valid correlation functions
                validCorrFunSlices = sum(isnan(corrFun(:,:,iBand1,iBand2)),1)==0;
                if sum(validCorrFunSlices)>2
                    [meanCC,CI] = correlationBootstrap(corrFun(:,validCorrFunSlices,iBand1,iBand2),...
                        bounds(1,validCorrFunSlices,iBand1,iBand2),p.nBoot,p.alpha);
                    bootstrapCorrFun(:,iBand1,iBand2)=meanCC;
                    bootstrapBounds(:,:,iBand1,iBand2)=CI;
                end   
                
            end
            if ishandle(wtBar), waitbar(i1/numel(bands1),wtBar); end
        end
        
        save(outFilePaths{iInput1,iInput2},'corrFun','bounds','lags',...
        'bootstrapCorrFun','bootstrapBounds');
    end
end

disp('Finished calculating correlation...')
if ishandle(wtBar)&&isempty(ip.Results.waitbar), close(wtBar); end

end