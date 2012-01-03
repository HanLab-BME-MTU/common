function calculateMovieCoherence(movieObject,varargin)
% calculateMovieCorrelation calculate the spectral density of temporal maps
%
% SYNOPSIS calculateMovieCoherence(movieObject,paramsIn)
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
% Sebastien Besson, Sep 2011
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
iProc = movieObject.getProcessIndex('CoherenceCalculationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieObject.processes_)+1;
    movieObject.addProcess(CoherenceCalculationProcess(movieObject,...
        movieObject.outputDirectory_));                                                                                                 
end

cohereProc = movieObject.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(cohereProc,paramsIn);


stack=dbstack;
if ~any(strcmp('Process.run',{stack(:).name}));
    cohereProc.run();
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
    nMovies = numel(p.MovieIndex);
    movieCohereProc=cell(nMovies,1);
    movieInput = cell(nMovies,1);
    
    for i =1:nMovies;
        iMovie = p.MovieIndex(i);
        fprintf(1,'Calculating coherence for movie %g/%g\n',i,nMovies);
        
        % Delegate  calculation for each movie of the list
        movieParams.SliceIndex=p.SliceIndex{iMovie};
        movieData = movieObject.movies_{iMovie};
        iProc = movieData.getProcessIndex('CoherenceCalculationProcess',1,0);
        if isempty(iProc)
            movieData.addProcess(CoherenceCalculationProcess(movieData,...
                movieData.outputDirectory_));
        end
        movieCohereProc{i} = movieData.processes_{end};
        parseProcessParams(movieCohereProc{i},movieParams);
        movieCohereProc{i}.run(wtBarArgs{:});
        
        movieInput{i}=movieCohereProc{i}.getInput;
        
    end  
    
    % Check input is the same for all movies
    assert(all(cellfun(@(x) isequal(x,movieInput{1}),movieInput)));
    input=movieInput{1};
    nInput= numel(movieInput{1});
    
    % Check number of frames per movie and determine lag limits
    nFrames =  cellfun(@(x) x.nFrames_,movieObject.movies_);
    nFrames=max(nFrames);
    timeInterval=unique(cellfun(@(x) x.timeInterval_,movieObject.movies_));
    assert(numel(timeInterval)==1);
    
    % Load input
    fprintf(1,'Calculating coherence for movie list\n');
    inFilePaths = cell(nInput,nMovies);
    data = cell(nInput,1);
    range = cell(nInput,1);
    disp('Using preprocessed signal from:');
    for i =1:nMovies;
        iMovie = p.MovieIndex(i);
        [inFilePaths(:,i),localdata,localrange] =getTSInput(movieObject.movies_{iMovie},movieInput{i});
        for iInput=1:nInput
            for iBand=1:min(numel(data{iInput}),numel(localdata{iInput}))
                data{iInput}{iBand} =vertcat(data{iInput}{iBand},localdata{iInput}{iBand});
                range{iInput}{iBand} =vertcat(range{iInput}{iBand},localrange{iInput}{iBand});
            end
            for iBand=numel(data{iInput})+1:max(numel(data{iInput}),numel(localdata{iInput}))
                data{iInput}{iBand} =localdata{iInput}{iBand};
                range{iInput}{iBand} =localrange{iInput}{iBand};
            end
        end
    end
    cohereProc.setInFilePaths(inFilePaths);
    p.SliceIndex = vertcat(p.SliceIndex{:});
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
    
    input = cohereProc.getInput;
    nInput=numel(input);
    disp('Using preprocessed signal from:');
    [inFilePaths,data,range] = getTSInput(movieData,input);
    cohereProc.setInFilePaths(inFilePaths);
end

% Set up output files
outFilePaths=cell(nInput,nInput);
for i=1:nInput
    for j=1:i-1
        outFilePaths{i,j} = [p.OutputDirectory filesep 'coherence' ...
            input(i).name '_' input(j).name '.mat'];
    end
    outFilePaths{i,i} = [p.OutputDirectory filesep 'powerSpectrum' ...
        input(i).name '.mat'];
end
mkClrDir(p.OutputDirectory);
cohereProc.setOutFilePaths(outFilePaths);
disp('Results will be saved under:')
disp(p.OutputDirectory);

%% --------------- Coherence calculation ---------------%%% 

%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nfft = 2^nextpow2(nFrames); % cf pwelch
nFreqMax=nfft/2+1;
fs =1/timeInterval;
f= fs/2*linspace(0,1,nfft/2 +1); %#ok<NASGU>
nBands =cellfun(@numel,data);

padTS = @(x) padarray(x,nFrames-length(x),0,'post');

logMsg = @(i,j) ['Please wait, calculating ' input(i).name '/'...
    input(j).name ' coherence'];

% Calculate spectral density coherence
for iInput1=1:nInput
    for iInput2=1:iInput1-1
        disp(logMsg(iInput1,iInput2));
        
        % Initialize cross-correlation function and bounds
        avgCoh = NaN(nFreqMax,nBands(iInput1),nBands(iInput2));
        cohCI = NaN(2,nFreqMax,nBands(iInput1),nBands(iInput2));
        
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput1,iInput2)); end
        
        % Loop over bands and window slices
        bands1=p.BandMin:min(nBands(iInput1),p.BandMax);
        for i1=1:numel(bands1)
            iBand1 = bands1(i1);
            for iBand2=p.BandMin:min(nBands(iInput2),p.BandMax)
                
                % Find valid range and test minimum number of timepoints
                nTimepoints = cellfun(@(x,y) length(intersect(x,y)),range{iInput2}{iBand2},...
                    range{iInput1}{iBand1});
                validSlices = nTimepoints>=minP & p.SliceIndex;
                
                if sum(validSlices)>0
                    % Concatenate time-series
                    paddedTS1 = cellfun(padTS,data{iInput1}{iBand1}(validSlices),'Unif',false);
                    paddedTS1 = cat(2,paddedTS1{:});
                    paddedTS2 = cellfun(padTS,data{iInput2}{iBand2}(validSlices),'Unif',false);
                    paddedTS2 = cat(2,paddedTS2{:});
                    
                    % Bootstrap coherence
                    [c,cI]=coherenceBootstrap(paddedTS1,paddedTS2,p.nWin,...
                        p.window,p.noLap,fs,'alpha',p.alpha,'nBoot',p.nBoot);
                    avgCoh(:,iBand1,iBand2)=c;
                    cohCI(:,:,iBand1,iBand2)=cI;
                end
            end
            if ishandle(wtBar), waitbar(i1/numel(bands1),wtBar); end
        end
        
        save(outFilePaths{iInput1,iInput2},'avgCoh','cohCI','f');
    end
end

disp('Finished calculating coherence...')
if ishandle(wtBar)&&isempty(ip.Results.waitbar), close(wtBar); end

end

function [paths,data,range] = getTSInput(movieData,input)

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
    error(['Each time series must have been preprocessed !' ...
        'Please apply pre-processing to all time series before '...
        'running coherence calculation!'])
end

disp(signalPreproc.funParams_.OutputDirectory);

% Load input
paths = cell(nInput,1);
data = cell(nInput,1);
range = cell(nInput,1);
for iInput=1:nInput
    paths{iInput,1} = signalPreproc.outFilePaths_{1,preprocIndex(iInput)};
    [data{iInput},range{iInput}] = signalPreproc.loadChannelOutput(preprocIndex(iInput));
end
end