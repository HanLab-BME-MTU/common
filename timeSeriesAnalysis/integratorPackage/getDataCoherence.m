function outFilePaths = getDataCoherence(input,p,p2,varargin)
% Check input
ip=inputParser;
ip.addRequired('input',@isstruct);
ip.addRequired('p',@isstruct);
ip.addRequired('p2',@isstruct);
ip.addParamValue('waitbar',-1,@ishandle);
ip.parse(input,p,p2,varargin{:})

% Retrieve waitbar or create one if necessary
if ~isempty(ip.Results.waitbar)
    wtBar=ip.Results.waitbar;
elseif feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...');
else
    wtBar=-1;
end

% Set up output files
nInput=numel(input);
outFilePaths=cell(nInput,nInput);
for i=1:nInput
    for j=1:i-1
        outFilePaths{i,j} = [p2.outputDir filesep 'coherence' ...
            input(i).name '_' input(j).name '.mat'];
    end
    outFilePaths{i,i} = [p2.outputDir filesep 'powerSpectrum' ...
        input(i).name '.mat'];
end
disp('Starting calculating coherence...')
disp('Saving results under:');
disp(p2.outputDir);

%% Coherence calculation
%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nfft = 2^nextpow2(p.nFrames); % cf pwelch
nFreqMax=nfft/2+1;
fs =1/p.timeInterval;
f= fs/2*linspace(0,1,nfft/2 +1); %#ok<NASGU>
data={input.data};
range={input.range};
nBands =cellfun(@numel,data);

padTS = @(x) padarray(x,p.nFrames-length(x),0,'post');

logMsg = @(i,j) ['Calculating ' input(i).name '/' input(j).name ' coherence'];

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
                    [c,cI]=coherenceBootstrap(paddedTS1,paddedTS2,p2.nWin,...
                        p2.window,p2.noLap,fs,'alpha',p2.alpha,'nBoot',p2.nBoot);
                    avgCoh(:,iBand1,iBand2)=c;
                    cohCI(:,:,iBand1,iBand2)=cI;
                end
            end
            if ishandle(wtBar), waitbar(i1/numel(bands1),wtBar); end
        end
        
        save(outFilePaths{iInput1,iInput2},'f','avgCoh','cohCI');
    end
end

disp('Finished calculating coherence...')
