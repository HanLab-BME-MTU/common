function getMovieCorrelation(movieObject,data,range,varargin)

% Check input
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

%%
if isa(movieObject,'MovieList')
    nFrames = unique(cellfun(@(x) x.nFrames_,movieObject.movies_));
    timeInterval = unique(cellfun(@(x) x.timeInterval_,movieObject.movies_));
    assert(isscalar(nFrames) && isscalar(timeInterval));
else
    nFrames=movieObject.nFrames_;
    timeInterval =movieObject.timeInterval_;
end
assert(~isempty(nFrames) && ~isempty(timeInterval));

input = signalProc.getInput;
nInput=numel(input);

% Test the presence and output validity of the signal preprocessing
iSignalPreproc =movieObject.getProcessIndex('SignalPreprocessingProcess',1,1);
if isempty(iSignalPreproc)
    error([SignalPreprocessingProcess.getName ' has not yet been performed'...
        'on this movie! Please run first!!']);
end

% Check that there is a valid output
signalPreproc = movieObject.processes_{iSignalPreproc};
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
        'running correlation calculation!'])
end



%% --------------- Correlation calculation ---------------%%% 

%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nInput =numel(data);
nLagsMax =round(movieObject.nFrames_/4);
nBands =cellfun(@numel,data);
nSlices = numel(data{1}{1});

logMsg = @(i) ['Calculating ' input(i).name ' autocorrelation'];

% Calculate autocorrelation
lags =(0:nLagsMax)'*process.owner_.timeInterval; %#ok<NASGU>
for iInput=1:nInput
    disp(logMsg(iInput));
    
    % Initialize autocorrelation function and bounds
    acf = NaN(nLagsMax+1,nSlices,nBands(iInput));
    acfBounds = NaN(2,nSlices,nBands(iInput));
    bootstrapAcf=NaN(nLagsMax+1,nBands(iInput));
    bootstrapAcfBounds=NaN(2,nLagsMax+1,nBands(iInput));
  
    pacf = NaN(nLagsMax+1,nSlices,nBands(iInput));
    pacfBounds = NaN(2,nSlices,nBands(iInput));
    bootstrapPacf=NaN(nLagsMax+1,nBands(iInput));
    bootstrapPacfBounds=NaN(2,nLagsMax+1,nBands(iInput));
    
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
    
    for iBand=p.BandMin:min(nBands(iInput),p.BandMax);
        
        % Get number of timepoints and prune out slices
        nTimepoints=cellfun(@length,data{iInput}{iBand});
        validSlices =nTimepoints >=minP & p.SliceIndex;
        
        % Calculate raw auto-correlation
        for iSlice=find(validSlices)'
            nLags = round(length(data{iInput}{iBand}{iSlice})/4);
            [acf(1:nLags+1,iSlice,iBand),~,acfBounds(:,iSlice,iBand)] = ...
                autocorr(data{iInput}{iBand}{iSlice},nLags);
            
            [pacf(1:nLags+1,iSlice,iBand),~,pacfBounds(:,iSlice,iBand)] = ...
                parcorr(data{iInput}{iBand}{iSlice},nLags);
        end
        
        % Bootstrap valid autocorrelation functions
        validAcfSlices = sum(isnan(acf(:,:,iBand)),1)==0;
        if sum(validAcfSlices)>2
            [meanCC,CI] = correlationBootstrap(acf(:,validAcfSlices,iBand),...
                acfBounds(1,validAcfSlices,iBand),p.nBoot,p.alpha);
            bootstrapAcf(:,iBand)=meanCC;
            bootstrapAcfBounds(:,:,iBand)=CI;
        end
        
        % Bootstrap valid partial autocorrelation functions
        validPacfSlices = sum(isnan(pacf(:,:,iBand)),1)==0;
        if sum(validPacfSlices)>2  
            [meanCC,CI] = correlationBootstrap(pacf(:,validPacfSlices,iBand),...
                pacfBounds(1,validPacfSlices,iBand),p.nBoot,p.alpha);
            bootstrapPacf(:,iBand)=meanCC;
            bootstrapPacfBounds(:,:,iBand)=CI;
        end
        
        if ishandle(wtBar), waitbar(iBand/nBands(iInput),wtBar); end
    end
    
    save(outFilePaths{iInput,iInput},'lags','acf','acfBounds',...
        'bootstrapAcf','bootstrapAcfBounds','pacf','pacfBounds',...
        'bootstrapPacf','bootstrapPacfBounds','-append');  
end

logMsg = @(i,j) ['Calculating ' input(i).name '/' input(j).name ' cross-correlation'];

% Calculate cross-correlation
lags =(-nLagsMax:nLagsMax)'*movieData.timeInterval_; %#ok<NASGU>
for iInput1=1:nInput
    for iInput2=1:iInput1-1
        disp(logMsg(iInput1,iInput2));
        
        % Initialize cross-correlation function and bounds
        ccf = NaN(2*nLagsMax+1,nSlices,nBands(iInput1),nBands(iInput2));
        ccfBounds  = NaN(2,nSlices,nBands(iInput1),nBands(iInput2));
        bootstrapCcf=NaN(2*nLagsMax+1,nBands(iInput1),nBands(iInput2));
        bootstrapCcfBounds=NaN(2,2*nLagsMax+1,nBands(iInput1),nBands(iInput2));
        
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
                    [ccf(1:2*nLags+1,iSlice,iBand1,iBand2),~,ccfBounds(:,iSlice,iBand1,iBand2)] =...
                        crosscorr(data{iInput1}{iBand1}{iSlice}(range1),data{iInput2}{iBand2}{iSlice}(range2),nLags);
                end
                
                % Bootstrap valid correlation functions
                validCcfSlices = sum(isnan(ccf(:,:,iBand1,iBand2)),1)==0;
                if sum(validCcfSlices)>2
                    [meanCC,CI] = correlationBootstrap(ccf(:,validCcfSlices,iBand1,iBand2),...
                        ccfBounds(1,validCcfSlices,iBand1,iBand2),p.nBoot,p.alpha);
                    bootstrapCcf(:,iBand1,iBand2)=meanCC;
                    bootstrapCcfBounds(:,:,iBand1,iBand2)=CI;
                end   
                
            end
            if ishandle(wtBar), waitbar(i1/numel(bands1),wtBar); end
        end
        
        save(outFilePaths{iInput1,iInput2},'lags','ccf','ccfBounds',...
        'bootstrapCcf','bootstrapCcfBounds','-append');
    end
end