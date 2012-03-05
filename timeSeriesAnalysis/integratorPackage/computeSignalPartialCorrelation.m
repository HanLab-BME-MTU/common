function outFilePaths = computeSignalPartialCorrelation(signal,params,varargin)

% Check input
ip=inputParser;
ip.addRequired('signal',@isstruct);
ip.addRequired('params',@isstruct);
ip.addParamValue('waitbar',-1,@ishandle);
ip.parse(signal,params,varargin{:})

% Retrieve waitbar or create one if necessary
if ~isempty(ip.Results.waitbar)
    wtBar=ip.Results.waitbar;
elseif feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...');
else
    wtBar=-1;
end

% Set up output files
nSignal=numel(signal);
outFilePaths=cell(nSignal,nSignal);
for i=1:nSignal 
    outFilePaths{i,i} = [params.OutputDirectory filesep 'partialautocorrelation_' ...
        signal(i).name '.mat'];
end
disp('Starting calculating partial correlation...')
disp('Saving results under:');
disp(params.OutputDirectory);

%% Partial correlation calculation
%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nSignal =numel(signal);
nLagsMax =round(params.nFrames/4);
data={signal.data};
nBands =cellfun(@numel,data);
nSlices = numel(data{1}{1});


logMsg = @(i) ['Calculating ' signal(i).name ' partial autocorrelation'];

% Calculate autocorrelation
lags =(0:nLagsMax)'*params.timeInterval; %#ok<NASGU>
for iSignal=1:nSignal
    disp(logMsg(iSignal));
    
    pacf = NaN(nLagsMax+1,nSlices,nBands(iSignal));
    pacfBounds = NaN(2,nSlices,nBands(iSignal));
    bootstrapPacf=NaN(nLagsMax+1,nBands(iSignal));
    bootstrapPacfBounds=NaN(2,nLagsMax+1,nBands(iSignal));
    
    
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iSignal)); end
    
    for iBand=1:nBands(iSignal)
        % Get number of timepoints and prune out slices
        nTimepoints=cellfun(@length,data{iSignal}{iBand});
        validSlices =nTimepoints >=minP;
        
        % Calculate raw auto-correlation
        for iSlice=find(validSlices)'
            nLags = round(length(data{iSignal}{iBand}{iSlice})/4);
            [pacf(1:nLags+1,iSlice,iBand),~,pacfBounds(:,iSlice,iBand)] = ...
                parcorr(data{iSignal}{iBand}{iSlice},nLags);
        end
        
        
        % Bootstrap valid partial autocorrelation functions
        validPacfSlices = sum(isnan(pacf(:,:,iBand)),1)==0;
        if sum(validPacfSlices)>2  
            [meanCC,CI] = correlationBootstrap(pacf(:,validPacfSlices,iBand),...
                pacfBounds(1,validPacfSlices,iBand),params.nBoot,params.alpha);
            bootstrapPacf(:,iBand)=meanCC;
            bootstrapPacfBounds(:,:,iBand)=CI;
        end
        
        if ishandle(wtBar), waitbar(iBand/nBands(iSignal),wtBar); end
    end
    
    save(outFilePaths{iSignal,iSignal},'lags','pacf','pacfBounds',...
        'bootstrapPacf','bootstrapPacfBounds');  
end


disp('Finished calculating partial correlation...')

