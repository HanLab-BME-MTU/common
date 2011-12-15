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
    
end

assert(all(cellfun(@(x) isequal(x,movieInput{1}),movieInput)));

input=movieInput{1};
nCorrFun=nInput*(nInput+1)/2;
nFrames =  cellfun(@(x) x.nFrames_,movieList.movies_);
inFilePaths = cell(nCorrFun,nMovies);

rawCorrFun = cell(nCorrFun,1);
rawBounds = cell(nMovies,1);
corrFun=cell(nCorrFun,1);
bounds=cell(nCorrFun,1);

nAcfLags=max(nFrames)/4+1;
nCcfLags=max(nFrames)/2+1;
bandRange= p.BandMin:p.BandMax;
for iInput=1:nInput
    for jInput=1:iInput
        iCorrFun = iInput*(iInput-1)/2+ jInput;
        for i=1:nMovies
            iMovie=p.MovieIndex(i);
            inFilePaths{iCorrFun,i}=movieCorrProc(i).outFilePaths_{iInput,jInput};
            [rawCorrFun{iCorrFun,i} rawBounds{iCorrFun,i}] =...
                movieCorrProc(i).loadChannelOutput(iInput,jInput,'output',{'corrFun','bounds'});
            
            % Slice window slices
            rawCorrFun{iCorrFun,i}=rawCorrFun{iCorrFun,i}(:,find(p.SliceIndex{iMovie}),:,:);
            rawBounds{iCorrFun,i}=rawBounds{iCorrFun,i}(:,find(p.SliceIndex{iMovie}),:,:);
            
            % Slice bands
            if size(rawCorrFun{iCorrFun,i},3)>1
                rawCorrFun{iCorrFun,i}=rawCorrFun{iCorrFun,i}(:,:,bandRange,:);
                rawBounds{iCorrFun,i}=rawBounds{iCorrFun,i}(:,:,bandRange,:);
            end
            if size(rawCorrFun{iCorrFun,i},4)>1
                rawCorrFun{iCorrFun,i}=rawCorrFun{iCorrFun,i}(:,:,:,bandRange);
                rawBounds{iCorrFun,i}=rawBounds{iCorrFun,i}(:,:,:,bandRange);
            end
            
            % Append
            if iInput==jInput
                rawCorrFun{iCorrFun,i}(end+1:nAcfLags,:,:,:)=NaN;
            else
                rawCorrFun{iCorrFun,i}(end+1:nCcfLags,:,:,:)=NaN;
            end
        end
        corrFun{iCorrFun}=horzcat([rawCorrFun{iCorrFun,:}]);
        bounds{iCorrFun}=horzcat([rawBounds{iCorrFun,:}]);
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
mkClrDir(p.OutputDirectory);
corrProc.setOutFilePaths(outFilePaths);

%% --------------- Correlation calculation ---------------%%%
disp('Starting bootstraping correlation...')

acfLogMsg = @(i) ['Please wait, bootstrapping ' input(i).name ' autocorrelation'];
ccfLogMsg = @(i,j) ['Please wait, bootstrapping ' input(i).name '/'...
    input(j).name ' cross-correlation'];


% Calculate autocorrelation
for iInput1=1:nInput
    for iInput2=1:iInput1
        iCorrFun = iInput1*(iInput1-1)/2+ iInput2;
        if iInput1==iInput2,
            disp(acfLogMsg(iInput1));
        else
            disp(ccfLogMsg(iInput1,iInput2));
        end
        
        
        % Initialize autocorrelation function and bounds
        nBandMax1 =min(p.BandMax,size(corrFun{iCorrFun},3));
        nBandMax2 =min(p.BandMax,size(corrFun{iCorrFun},4));
        if iInput1==iInput2,
            lags=(0:nAcfLags-1)*movieList.movies_{1}.timeInterval_;
        else
            lags=(-(nCcfLags-1)/2:(nCcfLags-1)/2)*movieList.movies_{1}.timeInterval_;
        end
        bootstrapCorrFun=NaN(numel(lags),nBandMax1,nBandMax2);
        bootstrapBounds=NaN(2,numel(lags),nBandMax1,nBandMax2);
        
        if ishandle(wtBar),
            if iInput1==iInput2,
                waitbar(0,wtBar,acfLogMsg(iInput1));
            else
                waitbar(0,wtBar,ccfLogMsg(iInput1,iInput2));
            end
        end
        
        % Loop over bands and window slices
        for iBand1=p.BandMin:nBandMax1
            for iBand2=p.BandMin:nBandMax2
                
                % Bootstrap valid correlation functions
                validSlices = sum(isnan(corrFun{iCorrFun}(:,:,iBand1,iBand2)),1)==0;
                if sum(validSlices)>2
                    [meanCC,CI] = correlationBootstrap(corrFun{iCorrFun}(:,validSlices,iBand1,iBand2),...
                        bounds{iCorrFun}(1,validSlices,iBand1,iBand2),p.nBoot,p.alpha);
                    bootstrapCorrFun(:,iBand1,iBand2)=meanCC;
                    bootstrapBounds(:,:,iBand1,iBand2)=CI;
                end
                
            end
            
            if ishandle(wtBar), waitbar(iBand1/nBandMax1,wtBar); end
        end
        
        save(outFilePaths{iInput1,iInput2},'lags','bootstrapCorrFun','bootstrapBounds');
    end
end

disp('Finished bootstraping correlation...')
if ishandle(wtBar) && isempty(ip.Results.waitbar), close(wtBar); end

end