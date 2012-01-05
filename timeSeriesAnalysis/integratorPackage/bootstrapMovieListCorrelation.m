function bootstrapMovieListCorrelation(movieList,varargin)
% bootstrapMovieListCorrelation calculate the autocorrelation and cross-correlation
% between the protrusion and activity maps
%
% SYNOPSIS calculateMovieCorrelation(movieList,paramsIn)
%
% INPUT
%   movieList - A MovieList object
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Marco Vilela, Sep 2011
% Sebastien Besson, Dec 2011
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieList', @(x) isa(x,'MovieList'));
ip.addOptional('paramsIn',[], @isstruct);
ip.addParamValue('waitbar',[], @ishandle);
ip.parse(movieList,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes
iProc = movieList.getProcessIndex('CorrelationCalculationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieObject.processes_)+1;
    movieList.addProcess(CorrelationCalculationProcess(movieObject,...
        movieObject.outputDirectory_));
end

corrProc = movieList.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(corrProc,paramsIn);


%% --------------- Initialization ---------------%%
disp('Starting bootstrapping movie list correlation functions...')
[~,movieName]=fileparts(movieList.getPath);
if ~isempty(ip.Results.waitbar)
    wtBar=ip.Results.waitbar;
    waitbar(0,wtBar,'Initializing...','Name',movieName);
elseif feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',movieName);
else
    wtBar=-1;
end



% Load input
nMovies =  numel(p.MovieIndex);
movieInput= cell(nMovies,1);
movieCorrProc(nMovies,1)=CorrelationCalculationProcess;
disp('Using correlation functions from:');
for i=1:nMovies
    iMovie=p.MovieIndex(i);
    iMovieCorrProc =movieList.movies_{iMovie}.getProcessIndex('CorrelationCalculationProcess',1,false);
    movieCorrProc(i) = movieList.movies_{iMovie}.processes_{iMovieCorrProc};
    
    movieInput{i}=movieCorrProc(i).getInput;
    nInput=numel(movieInput{i});
    
    for iInput=1:nInput
        for jInput=1:iInput-1
            if ~movieCorrProc(i).checkChannelOutput(iInput,jInput)
                error(['Missing correlation output !' ...
                    'Please apply pre-processing to all time series before '...
                    'running correlation calculatino!'])
            end
        end
    end
    
    disp(movieCorrProc(i).funParams_.OutputDirectory);
end

% Check that all input are constistent
assert(all(cellfun(@(x) isequal(x,movieInput{1}),movieInput)));
input=movieInput{1};
nCorrFun=nInput*(nInput+1)/2;

% Check number of frames per movie and determine lag limits
nFrames =  cellfun(@(x) x.nFrames_,movieList.movies_);
nAcfLags=max(nFrames)/4+1;
nCcfLags=max(nFrames)/2+1;

% Initialize input files and movie correlation functions
inFilePaths = cell(nCorrFun,nMovies);
bandRange= p.BandMin:p.BandMax;

% Get autocorrelation functions
acf=cell(nInput,1);
acfBounds=cell(nInput,1);
pacf=cell(nInput,1);
pacfBounds=cell(nInput,1);
for iInput=1:nInput
    iCorrFun = iInput*(iInput-1)/2+ iInput;
    rawAcf = cell(nMovies,1);
    rawAcfBounds = cell(nMovies,1);
    rawPacf = cell(nMovies,1);
    rawPacfBounds = cell(nMovies,1);
    for i=1:nMovies
        % Load raw correlation functions from the movie process
        iMovie=p.MovieIndex(i);
        inFilePaths{iCorrFun,i}=movieCorrProc(i).outFilePaths_{iInput,iInput};
        [rawAcf{i} rawAcfBounds{i} rawPacf{i} rawPacfBounds{i}] =...
            movieCorrProc(i).loadChannelOutput(iInput,iInput,...
            'output',{'acf','acfBounds','pacf','pacfBounds'});

        % Select only window slices given by p.SliceIndex and bands given
        % by bandRange
        rawAcf{i}=extractSlice(rawAcf{i},p.SliceIndex{iMovie},bandRange);
        rawAcfBounds{i}=extractSlice(rawAcfBounds{i},p.SliceIndex{iMovie},bandRange);
        rawPacf{i}=extractSlice(rawPacf{i},p.SliceIndex{iMovie},bandRange);
        rawPacfBounds{i}=extractSlice(rawPacfBounds{i},p.SliceIndex{iMovie},bandRange);
        
        % Append trailing NaNs if movies have different number of frames
        rawAcf{i}(end+1:nAcfLags,:,:,:)=NaN;
        rawPacf{i}(end+1:nAcfLags,:,:,:)=NaN;
    end
    % Concatenate bounds and correlation functions
    acf{iInput}=horzcat([rawAcf{:}]);
    acfBounds{iInput}=horzcat([rawAcfBounds{:}]);
    pacf{iInput}=horzcat([rawPacf{:}]);
    pacfBounds{iInput}=horzcat([rawPacfBounds{:}]);
end

% Get cross-correlation functions
ccf=cell(nInput,nInput);
ccfBounds=cell(nInput,nInput);
for iInput=1:nInput   
    for jInput=1:iInput-1
        iCorrFun = iInput*(iInput-1)/2+ jInput;
        rawCcf = cell(nMovies,1);
        rawCcfBounds = cell(nMovies,1);

        for i=1:nMovies
            % Load raw correlation functions from the movie process
            iMovie=p.MovieIndex(i);
            inFilePaths{iCorrFun,i}=movieCorrProc(i).outFilePaths_{iInput,iInput};
            [rawCcf{i} rawCcfBounds{i}] =...
                movieCorrProc(i).loadChannelOutput(iInput,jInput,'output',{'ccf','ccfBounds'});
            
            % Select only window slices given by p.SliceIndex and bands given
            % by bandRange
            rawCcf{i}=extractSlice(rawCcf{i},p.SliceIndex{iMovie},bandRange);
            rawCcfBounds{i}=extractSlice(rawCcfBounds{i},p.SliceIndex{iMovie},bandRange);
            
            % Append trailing NaNs if movies have different number of frames
            rawCcf{i}(end+1:nCcfLags,:,:,:)=NaN;
        end
        % Concatenate bounds and correlation functions
        ccf{iInput,jInput}=horzcat([rawCcf{:}]);
        ccfBounds{iInput,jInput}=horzcat([rawCcfBounds{:}]);
    end
end
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

% Calculate autocorrelation
logMsg = @(i) ['Bootstrapping ' input(i).name ' autocorrelation'];
lags=(0:nAcfLags-1)*movieList.movies_{1}.timeInterval_;
for iInput=1:nInput
    disp(logMsg(iInput));
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
    
    % Initialize autocorrelation function and bounds
    nBandMax =min(p.BandMax,size(acf{iInput},3));
    
    % Initialize output
    bootstrapAcf=NaN(numel(lags),nBandMax);
    bootstrapAcfBounds=NaN(2,numel(lags),nBandMax);
    bootstrapPacf=NaN(numel(lags),nBandMax);
    bootstrapPacfBounds=NaN(2,numel(lags),nBandMax);
    
    % Loop over bands and window slices
    for iBand=p.BandMin:nBandMax
            
        % Bootstrap valid autocorrelation functions
        validAcfSlices = sum(isnan(acf{iInput}(:,:,iBand)),1)==0;
        if sum(validAcfSlices)>2
            [meanCC,CI] = correlationBootstrap(acf{iInput}(:,validAcfSlices,iBand),...
                acfBounds{iInput}(1,validAcfSlices,iBand),p.nBoot,p.alpha);
            bootstrapAcf(:,iBand)=meanCC;
            bootstrapAcfBounds(:,:,iBand)=CI;
        end
        
        % Bootstrap valid partial autocorrelation functions
        validPacfSlices = sum(isnan(pacf{iInput}(:,:,iBand)),1)==0;
        if sum(validPacfSlices)>2
            [meanCC,CI] = correlationBootstrap(pacf{iInput}(:,validPacfSlices,iBand),...
                pacfBounds{iInput}(1,validPacfSlices,iBand),p.nBoot,p.alpha);
            bootstrapPacf(:,iBand)=meanCC;
            bootstrapPacfBounds(:,:,iBand)=CI;
        end
        
        if ishandle(wtBar), waitbar(iBand/nBandMax,wtBar); end
    end
    
    save(outFilePaths{iInput,iInput},'lags','acf','acfBounds',...
        'bootstrapAcf','bootstrapAcfBounds','pacf','pacfBounds',...
        'bootstrapPacf','bootstrapPacfBounds');  
end
    

% Calculate cross-correlation
logMsg = @(i,j) ['Bootstrapping ' input(i).name '/' input(j).name ' cross-correlation'];
lags=(-(nCcfLags-1)/2:(nCcfLags-1)/2)*movieList.movies_{1}.timeInterval_;
for iInput1=1:nInput
    for iInput2=1:iInput1-1
        disp(logMsg(iInput1,iInput2));
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput1,iInput2)); end
        
        % Initialize autocorrelation function and bounds
        nBandMax1 =min(p.BandMax,size(ccf{iInput1,iInput2},3));
        nBandMax2 =min(p.BandMax,size(ccf{iInput1,iInput2},4));
        
        % Initialize output
        bootstrapCcf=NaN(numel(lags),nBandMax1,nBandMax2);
        bootstrapCcfBounds=NaN(2,numel(lags),nBandMax1,nBandMax2);
        
        % Loop over bands and window slices
        for iBand1=p.BandMin:nBandMax1
            for iBand2=p.BandMin:nBandMax2
                
                % Bootstrap valid correlation functions
                validCcfSlices = sum(isnan(ccf{iInput1,iInput2}(:,:,iBand1,iBand2)),1)==0;
                if sum(validCcfSlices)>2
                    [meanCC,CI] = correlationBootstrap(ccf{iInput1,iInput2}(:,validCcfSlices,iBand1,iBand2),...
                        ccfBounds{iInput1,iInput2}(1,validCcfSlices,iBand1,iBand2),p.nBoot,p.alpha);
                    bootstrapCcf(:,iBand1,iBand2)=meanCC;
                    bootstrapCcfBounds(:,:,iBand1,iBand2)=CI;
                end
            end
            
            if ishandle(wtBar), waitbar(iBand1/nBandMax1,wtBar); end
        end
        save(outFilePaths{iInput1,iInput2},'lags','ccf','ccfBounds',...
            'bootstrapCcf','bootstrapCcfBounds');
    end
end

disp('Finished bootstraping correlation...')
if ishandle(wtBar) && isempty(ip.Results.waitbar), close(wtBar); end

end

function f = extractSlice(f,sliceIndex,bandRange)

f=f(:,sliceIndex,:,:);

% Select only band within the valid band range
if size(f,3)>1, f=f(:,:,bandRange,:);  end
if size(f,4)>1,f=f(:,:,:,bandRange); end

end