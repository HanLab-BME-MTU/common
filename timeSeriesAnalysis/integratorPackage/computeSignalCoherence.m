function outFilePaths = computeSignalCoherence(signal,params,varargin)
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
        outFilePaths{i,j} = [params.OutputDirectory filesep 'coherence' ...
            signal(i).name '_' signal(j).name '.mat'];
    end
    outFilePaths{i,i} = [params.OutputDirectory filesep 'powerSpectrum' ...
        signal(i).name '.mat'];
end
disp('Starting calculating coherence...')
disp('Saving results under:');
disp(params.OutputDirectory);

%% Coherence calculation
%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nfft = 2^nextpow2(params.nFrames); % cf pwelch
nFreqMax=nfft/2+1;
fs =1/params.timeInterval;
f= fs/2*linspace(0,1,nfft/2 +1); %#ok<NASGU>
data={signal.data};
range={signal.range};
nBands =cellfun(@numel,data);

padTS = @(x) padarray(x,params.nFrames-length(x),0,'post');

logMsg = @(i,j) ['Calculating ' signal(i).name '/' signal(j).name ' coherence'];

% Calculate spectral density coherence
for iSignal1=1:nSignal
    for iSignal2=1:iSignal1-1
        disp(logMsg(iSignal1,iSignal2));
        
        % Initialize cross-correlation function and bounds
        avgCoh = NaN(nFreqMax,nBands(iSignal1),nBands(iSignal2));
        cohCI = NaN(2,nFreqMax,nBands(iSignal1),nBands(iSignal2));
        
        if ishandle(wtBar), waitbar(0,wtBar,logMsg(iSignal1,iSignal2)); end
        
        % Loop over bands and window slices
        for iBand1=1:nBands(iSignal1)
            for iBand2=1:nBands(iSignal2)
                
                % Find valid range and test minimum number of timepoints
                nTimepoints = cellfun(@(x,y) length(intersect(x,y)),range{iSignal2}{iBand2},...
                    range{iSignal1}{iBand1});
                validSlices = nTimepoints>=minP;
                
                if sum(validSlices)>0
                    % Concatenate time-series
                    paddedTS1 = cellfun(padTS,data{iSignal1}{iBand1}(validSlices),'Unif',false);
                    paddedTS1 = cat(2,paddedTS1{:});
                    paddedTS2 = cellfun(padTS,data{iSignal2}{iBand2}(validSlices),'Unif',false);
                    paddedTS2 = cat(2,paddedTS2{:});
                    
                    % Bootstrap coherence
                    [c,cI]=coherenceBootstrap(paddedTS1,paddedTS2,params.nWin,...
                        params.window,params.noLap,fs,'alpha',params.alpha,'nBoot',params.nBoot);
                    avgCoh(:,iBand1,iBand2)=c;
                    cohCI(:,:,iBand1,iBand2)=cI;
                end
            end
            if ishandle(wtBar), waitbar(iBand1/nBands(iSignal1),wtBar); end
        end
        
        save(outFilePaths{iSignal1,iSignal2},'f','avgCoh','cohCI');
    end
end

disp('Finished calculating coherence...')
