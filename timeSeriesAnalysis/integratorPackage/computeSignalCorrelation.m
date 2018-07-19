function outFilePaths = computeSignalCorrelation(signal,params,varargin)

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
    for j=1:i-1
        outFilePaths{i,j} = [params.OutputDirectory filesep 'crosscorrelation_' ...
            signal(i).name '_' signal(j).name '.mat'];
    end
    outFilePaths{i,i} = [params.OutputDirectory filesep 'autocorrelation_' ...
        signal(i).name '.mat'];
end
disp('Starting calculating correlation...')
disp('Saving results under:');
disp(params.OutputDirectory);

%% Correlation calculation 
%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nSignal =numel(signal);
nLagsMax =round(params.nFrames/4);
data={signal.data};
range={signal.range};
nBands =cellfun(@numel,data);
nSlices = numel(data{1}{1});

logMsg = @(i) ['Calculating ' signal(i).name ' autocorrelation'];

% Calculate autocorrelation
lags =(0:nLagsMax)'*params.timeInterval; %#ok<NASGU>
for iSignal=1:nSignal
    disp(logMsg(iSignal));
    
    % Initialize autocorrelation function and bounds
    acf = NaN(nLagsMax+1,nSlices,nBands(iSignal));
    acfBounds = NaN(2,nSlices,nBands(iSignal));
    bootstrapAcf=NaN(nLagsMax+1,nBands(iSignal));
    bootstrapAcfBounds=NaN(2,nLagsMax+1,nBands(iSignal));
    
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iSignal)); end
    
    for iBand=1:nBands(iSignal)
        % Get number of timepoints and prune out slices
        nTimepoints=cellfun(@length,data{iSignal}{iBand});
        validSlices =nTimepoints >=minP;
        
        % Calculate raw auto-correlation
        for iSlice=find(validSlices)'
            nLags = round(length(data{iSignal}{iBand}{iSlice})/4);
            [acf(1:nLags+1,iSlice,iBand),~,acfBounds(:,iSlice,iBand)] = ...
                autocorr(data{iSignal}{iBand}{iSlice},nLags);
        end
        
        % Bootstrap valid autocorrelation functions
        validAcfSlices = sum(isnan(acf(:,:,iBand)),1)==0;
        if sum(validAcfSlices)>2
            [meanCC,CI] = correlationBootstrap(acf(:,validAcfSlices,iBand),...Movie
                acfBounds(1,validAcfSlices,iBand),params.nBoot,params.alpha);
            bootstrapAcf(:,iBand)=meanCC;
            bootstrapAcfBounds(:,:,iBand)=CI;
        end
        
        if ishandle(wtBar), waitbar(iBand/nBands(iSignal),wtBar); end
    end
    
    save(outFilePaths{iSignal,iSignal},'lags','acf','acfBounds',...
        'bootstrapAcf','bootstrapAcfBounds');  
end

logMsg = @(i,j) ['Calculating ' signal(i).name '/' signal(j).name ' cross-correlation'];

% Calculate cross-correlation
lags =(-nLagsMax:nLagsMax)'*params.timeInterval; %#ok<NASGU>
for iSignal1=1:nSignal
    for iSignal2=1:iSignal1-1
        disp(logMsg(iSignal1,iSignal2));
        
        % Initialize cross-correlation function and bounds
        ccf = NaN(2*nLagsMax+1,nSlices,nBands(iSignal1),nBands(iSignal2));
        ccfBounds  = NaN(2,nSlices,nBands(iSignal1),nBands(iSignal2));
        bootstrapCcf=NaN(2*nLagsMax+1,nBands(iSignal1),nBands(iSignal2));
        bootstrapCcfBounds=NaN(2,2*nLagsMax+1,nBands(iSignal1),nBands(iSignal2));
        
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iSignal1,iSignal2)); end
        
        % Loop over bands and window slices
        for iBand1=1:nBands(iSignal1)
            for iBand2=1:nBands(iSignal2)
               
                % Find valid range and test minimum number of timepoints
                nTimepoints = cellfun(@(x,y) length(intersect(x,y)),range{iSignal2}{iBand2},...
                    range{iSignal1}{iBand1});
                validSlices = nTimepoints>=minP ;
                
                % Calculate raw cross-correlation
                for iSlice=find(validSlices)'
                    % Retrieve number of lags from range intersection
                    [~,range1,range2] = intersect(range{iSignal1}{iBand1}{iSlice},range{iSignal2}{iBand2}{iSlice});
                    nLags = round(length(range1)/4);
                    [ccf(1:2*nLags+1,iSlice,iBand1,iBand2),~,ccfBounds(:,iSlice,iBand1,iBand2)] =...
                        crosscorr(data{iSignal1}{iBand1}{iSlice}(range1),data{iSignal2}{iBand2}{iSlice}(range2),nLags);
                end
                
                % Bootstrap valid correlation functions
                validCcfSlices = sum(isnan(ccf(:,:,iBand1,iBand2)),1)==0;
                if sum(validCcfSlices)>2
                    [meanCC,CI] = correlationBootstrap(ccf(:,validCcfSlices,iBand1,iBand2),...
                        ccfBounds(1,validCcfSlices,iBand1,iBand2),params.nBoot,params.alpha);
                    bootstrapCcf(:,iBand1,iBand2)=meanCC;
                    bootstrapCcfBounds(:,:,iBand1,iBand2)=CI;
                end   
                
            end
            if ishandle(wtBar), waitbar(iBand1/nBands(iSignal1),wtBar); end
        end
        
        save(outFilePaths{iSignal1,iSignal2},'lags','ccf','ccfBounds',...
        'bootstrapCcf','bootstrapCcfBounds');
    end
end
disp('Finished calculating correlation...')

