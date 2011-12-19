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


%% --------------- Initialization ---------------%%
disp('Starting calculating spectral density...')
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
    for i =1:numel(p.MovieIndex);
        % Delegate correlation calculation for each movie of the list
        movieParams.SliceIndex=p.SliceIndex{p.MovieIndex(i)};
        movieData = movieObject.movies_{p.MovieIndex(i)};
        iProc = movieData.getProcessIndex('CoherenceCalculationProcess',1,0);
        if isempty(iProc)
            iProc = numel(movieData.processes_)+1;
            movieData.addProcess(CoherenceCalculationProcess(movieData,...
                movieData.outputDirectory_));
        end
        cohereProc = movieData.processes_{iProc};
        parseProcessParams(movieData.processes_{iProc},movieParams);
        cohereProc.run(wtBarArgs{:});
    end  
    
    % Calls the movie list correlation bootstrapping method
%     bootstrapMovieListCorrelation(movieObject,wtBarArgs{:});
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

input = cohereProc.getInput;
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
cohereProc.setInFilePaths(inFilePaths);

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

%% --------------- Correlation calculation ---------------%%% 
disp('Using preprocessed signal from:');
disp(signalPreproc.funParams_.OutputDirectory);
disp('Results will be saved under:')
disp(p.OutputDirectory);

%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nFreq = 256/2+1;
nBands =cellfun(@numel,data);
nSlices = numel(data{1}{1});

logMsg = @(i) ['Please wait, calculating ' input(i).name ' power spectral denisty'];

% Calculate autocorrelation
for iInput=1:nInput
    disp(logMsg(iInput));
    
    % Initialize spectral density
    P = NaN(nFreq,nSlices,nBands(iInput));P = NaN(nFreq,nSlices,nBands(iInput));
%     bootstrapP=NaN(nFreq,nBands(iInput));
        
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
    
    % Calculate valid band index
    validBands = find(~cellfun(@isempty,data{iInput}));
    for iBand=find(validBands<=p.BandMax & validBands>=p.BandMin )'
        
        % Calculate raw auto-correlation
        validSlices = ~cellfun(@isempty,data{iInput}{iBand});
        for iSlice=find(validSlices & p.SliceIndex)'
            P(:,iSlice,iBand) = pwelch(data{iInput}{iBand}{iSlice});
        end
        
%         % Bootstrap valid autocorrelation functions
%         validSlices = sum(isnan(corrFun(:,:,iBand)),1)==0;
%         if sum(validSlices)>2
%             [meanCC,CI] = correlationBootstrap(corrFun(:,validSlices,iBand),...
%                 bounds(1,validSlices,iBand),p.nBoot,p.alpha);
%             bootstrapCorrFun(:,iBand)=meanCC;
%             bootstrapBounds(:,:,iBand)=CI;
%         end
        
        if ishandle(wtBar), waitbar(iBand/nBands(iInput),wtBar); end
    end
       
    w =(0:pi/128:pi)'/movieData.timeInterval_;
    
    save(outFilePaths{iInput,iInput},'P','w');  
end

logMsg = @(i,j) ['Please wait, calculating ' input(i).name '/'...
    input(j).name ' coherence'];

% Calculate spectral density coherence
for iInput1=1:nInput
    for iInput2=1:iInput1-1
        disp(logMsg(iInput1,iInput2));
        
        % Initialize cross-correlation function and bounds
        P = NaN(nFreq,nSlices,nBands(iInput1),nBands(iInput2));
%         bootstrapCorrFun=NaN(2*nLagsMax+1,nBands(iInput1),nBands(iInput2));
        
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput1,iInput2)); end
        
        % Loop over bands and window slices
        for iBand1=p.BandMin:min(nBands(iInput1),p.BandMax)
            for iBand2=p.BandMin:min(nBands(iInput2),p.BandMax)
                
                % Calculate raw cross-correlation
                for iSlice=find(p.SliceIndex)'
                    % Find valid range and test minimum number of timepoints
                    [~,range1,range2] = intersect(range{iInput1}{iBand1}{iSlice},range{iInput2}{iBand2}{iSlice});
                    ccL               = length(range1);
                    if ccL >= minP
                        P(:,iSlice,iBand1,iBand2) =...
                            mscohere(data{iInput1}{iBand1}{iSlice}(range1),data{iInput2}{iBand2}{iSlice}(range2));
                    end
                end
                
                % Bootstrap valid correlation functions
%                 validSlices = sum(isnan(corrFun(:,:,iBand1,iBand2)),1)==0;
%                 if sum(validSlices)>2
%                     [meanCC,CI] = correlationBootstrap(corrFun(:,validSlices,iBand1,iBand2),...
%                         bounds(1,validSlices,iBand1,iBand2),p.nBoot,p.alpha);
%                     bootstrapCorrFun(:,iBand1,iBand2)=meanCC;
%                     bootstrapBounds(:,:,iBand1,iBand2)=CI;
%                 end   
%                 
            end
            if ishandle(wtBar), waitbar(iBand1/nBands(iInput1),wtBar); end
        end

        w =(-pi:pi/128:pi)'*movieData.timeInterval_;
        
        save(outFilePaths{iInput1,iInput2},'P','w');
    end
end

disp('Finished calculating spectral density...')
if ishandle(wtBar)&&isempty(ip.Results.waitbar), close(wtBar); end

end