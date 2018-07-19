function preprocessMovieSignal(movieObjectOrProcess,varargin)
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
%       ('BandMin' -> positive integer) Specifies the minimum band to
%       be used for pre-processing. Default 1.
%
%       ('BandMax' -> array of integers) Specified the last band to be
%       used for pre-processing the signal.
%
%       ('SliceIndex' -> logical array or cell array of logical array)  
%       If movieObject is a movie list, each element of the cell array 
%       specifies the SliceIndex of the slices to process for each movie. 
%       If movieObject is a movie, only the slices given by to true elements 
%       of the logical array will be processed.
%
%       ('kSigma'-> array of positive integer) For each input signal, the 
%       threshold (in sigma) to detect and remove time-series outliers. 
%       See detectOutliers.m for more info.
%
%       ('trendType'-> array of integer) For each input signal, a value 
%       specifying the type of detrending operation to apply. 
%       See removeMeanTrendNan for a list of acceptable values.
%
%       ('nBoot'-> integer) A value specifying the number of bootstrapped
%       samples to use when calculating the confidence interval of the
%       signal energy
%
%       ('alpha'-> scalar) A value specifying the alpha value to use when
%       calculating the confidence interval of the signal energy.

% Marco Vilela, Sep 2011
% Sebastien Besson, Nov 2011 (last modified Feb 2012)

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObjectOrProcess', @isProcessOrMovieObject);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieObjectOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get MovieData object and Process
[movieObject, signalPreProc] = getOwnerAndProcess(movieObjectOrProcess,'SignalPreprocessingProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(signalPreProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    [~,movieName]=fileparts(movieObject.getPath);
    wtBar = waitbar(0,'Initializing...','Name',movieName);
else
    wtBar=-1;
end

% If movie list, delegates signal preprocessing to individual movies
if isa(movieObject,'MovieList')
    movieParams=rmfield(p,{'MovieIndex','OutputDirectory'});
    nMovies =numel(p.MovieIndex);
    movieInput =cell(nMovies,1);
    movieSignalPreProc = cell(nMovies,1);
    for i =1:nMovies;
        % Get slice index 
        iMovie=p.MovieIndex(i);
        movieParams.SliceIndex=p.SliceIndex{iMovie};
        movieData = movieObject.getMovies{iMovie};
        fprintf(1,'Preprocessing signal for movie %g/%g\n',i,nMovies);
        
        [~,movieSignalPreProc{i}] = getOwnerAndProcess(movieData,'SignalPreprocessingProcess',true);
               
        % Parse parameters and run movie process
        parseProcessParams(movieSignalPreProc{i},movieParams);
        movieSignalPreProc{i}.run();
    end  
    
    % Check input is the same for all movies
    input=movieObject.getMovies{p.MovieIndex(1)}.getSampledOutput;
    nInput= numel(input);

    % Set input path
    paths=cellfun(@(x) x.outFilePaths_,movieSignalPreProc,'Unif',false);
    inFilePaths=vertcat(paths{:});
    signalPreProc.setInFilePaths(inFilePaths);
    
    % Load input for all  movies and concatenate them
    fprintf(1,'Preprocessing signal for movie list.\n');
    [input.energyData] = deal({});
    disp('Using preprocessed signal from:');
    for i =1:nMovies
        disp(inFilePaths{i});
        % Retrieve individual movie filtered output and concatenate all data
        for iInput=1:nInput
            localEnergyData = movieSignalPreProc{i}.loadOutput(iInput,'output','energyData');
            for iBand=1:min(numel(input(iInput).energyData),numel(localEnergyData))
                input(iInput).energyData{iBand} =horzcat(input(iInput).energyData{iBand},localEnergyData{iBand});
            end
            for iBand=numel(input(iInput).energyData)+1:max(numel(input(iInput).energyData),numel(localEnergyData))
                input(iInput).energyData{iBand} =localEnergyData{iBand};
            end
        end
    end
  
else  
    % Load sampled output
    input = movieObject.getSampledOutput;
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
end

% Set up output files
outFilePaths=cell(1,nInput);
for i=1:nInput
    outFilePaths{1,i} = [p.OutputDirectory filesep input(i).name '.mat'];
end
disp('Saving results under:');
disp(p.OutputDirectory);
mkClrDir(p.OutputDirectory);
signalPreProc.setOutFilePaths(outFilePaths);

% Create log messages
logMsg = @(iInput,iBand) ['Please wait, preprocessing signal of ' input(iInput).name];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
    
if isa(movieObject,'MovieList')
    %% Movie list energy calculation
    tic;
    % Initialize output
    nBands = arrayfun(@(x)numel(x.energyData),input);
    nBandsTot = sum(nBands);
    
    for iInput=1:nInput
        % Show log
        disp(logMsg(iInput));
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
        
        % Initialize output
        energy = NaN(nBands(iInput),3);
        
        % Calculate energy per input per band for all movies
        validBands=cellfun(@numel,input(iInput).energyData)>=2;
        for iBand=find(validBands);
            [energy(iBand,1), energy(iBand,2:3)] = ...
                getWindowsBandEnergy(input(iInput).energyData{iBand},p.nBoot,p.alpha);
            % Update waitbar
            if ishandle(wtBar),
                tj=toc;
                nj = sum(nBands(1:iInput-1))+ iBand;
                waitbar(nj/nBandsTot,wtBar,sprintf([logMsg(iInput) timeMsg(tj*nBandsTot/nj-tj)]));
            end
            
        end
        save(outFilePaths{1,iInput},'energy');
    end
else
    %% Movie signal preprocessing
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
    %minP     = 50; %The removeMeanTrendNaN code already checks this...
    
    % Create log messages
    tic;
    nBands=arrayfun(@(x) min(x,p.BandMax),nBands);
    nBandsTot = sum(nBands);
    for iInput=1:nInput
        % Show log
        disp(logMsg(iInput));
        if ishandle(wtBar), waitbar(sum(nBands(1:iInput-1))/nBandsTot ,wtBar,logMsg(iInput)); end
        
        % Initialize data and range ouptut
        rawData= cell(nBands(iInput),1);
        energyData = cell(nBands(iInput),1);
        energyRange = cell(nBands(iInput),1);
        data= cell(nBands(iInput),1);
        [data{:}]=deal(cell(nSlices,1));
        range = data;
        energy  = NaN(nBands(iInput),3);
        
        for iBand=p.BandMin:nBands(iInput)
            % Get iBand data and remove outliers
            rawData{iBand} =squeeze(inData{iInput}(:,iBand,:));
            rawData{iBand}(detectOutliers(rawData{iBand},p.kSigma(iInput))) = NaN;
            
            % Check percentage of NaN and remove linear trend for energy
            validSlices = p.SliceIndex;
            [energyData{iBand} ,energyRange{iBand},~,~,excVar] = ...
                removeMeanTrendNaN(rawData{iBand}(validSlices,:)','trendType',1);
            validSlices(excVar) = false;
            
            % Recompute signal detrending if using a different trend type
            if p.trendType(iInput)~=1
                [data{iBand}(validSlices) ,range{iBand}(validSlices)] = ...
                    removeMeanTrendNaN(rawData{iBand}(validSlices,:)','trendType',p.trendType(iInput));
            else
                data{iBand}(validSlices)=energyData{iBand};
                range{iBand}(validSlices)=energyRange{iBand};
            end
            
            % Calculate energy
            if sum(validSlices)>2, [energy(iBand,1), energy(iBand,2:3)] = getWindowsBandEnergy(energyData{iBand},p.nBoot,p.alpha); end
            
            % Update waitbar
            if ishandle(wtBar),
                tj=toc;
                nj = sum(nBands(1:iInput-1))+ iBand-p.BandMin+1;
                waitbar(nj/nBandsTot,wtBar,sprintf([logMsg(iInput) timeMsg(tj*nBandsTot/nj-tj)]));
            end
        end
        save(outFilePaths{1,iInput},'data','range','energyData','energy');
    end
    
    disp('Finished preprocessing signal...')
end

  
if ishandle(wtBar), close(wtBar); end
