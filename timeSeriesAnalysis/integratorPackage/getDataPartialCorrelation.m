function outFilePaths = getDataPartialCorrelation(input,p,p2,varargin)

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
   
    outFilePaths{i,i} = [p2.outputDir filesep 'partialautocorrelation_' ...
        input(i).name '.mat'];
end
disp('Starting calculating partial correlation...')
disp('Saving results under:');
disp(p2.outputDir);

%% Partial correlation calculation
%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
minP     = 50;

nInput =numel(input);
nLagsMax =round(p.nFrames/4);
data={input.data};
nBands =cellfun(@numel,data);
nSlices = numel(data{1}{1});


logMsg = @(i) ['Calculating ' input(i).name ' partial autocorrelation'];

% Calculate autocorrelation
lags =(0:nLagsMax)'*p.timeInterval; %#ok<NASGU>
for iInput=1:nInput
    disp(logMsg(iInput));
    
    pacf = NaN(nLagsMax+1,nSlices,nBands(iInput));
    pacfBounds = NaN(2,nSlices,nBands(iInput));
    bootstrapPacf=NaN(nLagsMax+1,nBands(iInput));
    bootstrapPacfBounds=NaN(2,nLagsMax+1,nBands(iInput));
    
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput)); end
    
    bandIndex = p.BandMin:min(nBands(iInput),p.BandMax);
    for i=1:numel(bandIndex)
        iBand=bandIndex(i);
        
        % Get number of timepoints and prune out slices
        nTimepoints=cellfun(@length,data{iInput}{iBand});
        validSlices =nTimepoints >=minP & p.SliceIndex;
        
        % Calculate raw auto-correlation
        for iSlice=find(validSlices)'
            nLags = round(length(data{iInput}{iBand}{iSlice})/4);
            [pacf(1:nLags+1,iSlice,iBand),~,pacfBounds(:,iSlice,iBand)] = ...
                parcorr(data{iInput}{iBand}{iSlice},nLags);
        end
        
        
        % Bootstrap valid partial autocorrelation functions
        validPacfSlices = sum(isnan(pacf(:,:,iBand)),1)==0;
        if sum(validPacfSlices)>2  
            [meanCC,CI] = correlationBootstrap(pacf(:,validPacfSlices,iBand),...
                pacfBounds(1,validPacfSlices,iBand),p2.nBoot,p2.alpha);
            bootstrapPacf(:,iBand)=meanCC;
            bootstrapPacfBounds(:,:,iBand)=CI;
        end
        
        if ishandle(wtBar), waitbar(iBand/numel(bandIndex),wtBar); end
    end
    
    save(outFilePaths{iInput,iInput},'lags','pacf','pacfBounds',...
        'bootstrapPacf','bootstrapPacfBounds');  
end

% logMsg = @(i,j) ['Calculating ' input(i).name '/' input(j).name ' partial cross-correlation'];
% 
% % Calculate cross-correlation
% lags =(-nLagsMax:nLagsMax)'*timeInterval; %#ok<NASGU>
% for iInput1=1:nInput
%     for iInput2=1:iInput1-1
%         disp(logMsg(iInput1,iInput2));
%         
%         % Initialize cross-correlation function and bounds
%         bootstrapCcf=NaN(2*nLagsMax+1,nBands(iInput1),nBands(iInput2));
%         bootstrapCcfBounds=NaN(2,2*nLagsMax+1,nBands(iInput1),nBands(iInput2));
%         
%         if ishandle(wtBar), waitbar(0,wtBar,logMsg(iInput1,iInput2)); end
%         
%         % Loop over bands and window slices
%         bands1=p.BandMin:min(nBands(iInput1),p.BandMax);
%         for i1=1:numel(bands1)
%             iBand1=bands1(i1);
%             for iBand2=p.BandMin:min(nBands(iInput2),p.BandMax)
%                
%                 % Find valid range and test minimum number of timepoints
%                 nTimepoints = cellfun(@(x,y) length(intersect(x,y)),range{iInput2}{iBand2},...
%                     range{iInput1}{iBand1});
%                 validSlices = nTimepoints>=minP & p.SliceIndex;
%                 
%                 % Calculate raw cross-correlation
%                 for iSlice=find(validSlices)'
%                     % Retrieve number of lags from range intersection
%                     [~,range1,range2] = intersect(range{iInput1}{iBand1}{iSlice},range{iInput2}{iBand2}{iSlice});
%                     nLags = round(length(range1)/4);
%                     [ccf(1:2*nLags+1,iSlice,iBand1,iBand2),~,ccfBounds(:,iSlice,iBand1,iBand2)] =...
%                         crosscorr(data{iInput1}{iBand1}{iSlice}(range1),data{iInput2}{iBand2}{iSlice}(range2),nLags);
%                 end
%                 
%                 % Bootstrap valid correlation functions
%                 validCcfSlices = sum(isnan(ccf(:,:,iBand1,iBand2)),1)==0;
%                 if sum(validCcfSlices)>2
%                     [meanCC,CI] = correlationBootstrap(ccf(:,validCcfSlices,iBand1,iBand2),...
%                         ccfBounds(1,validCcfSlices,iBand1,iBand2),p.nBoot,p.alpha);
%                     bootstrapCcf(:,iBand1,iBand2)=meanCC;
%                     bootstrapCcfBounds(:,:,iBand1,iBand2)=CI;
%                 end   
%                 
%             end
%             if ishandle(wtBar), waitbar(i1/numel(bands1),wtBar); end
%         end
%         
%         save(outFilePaths{iInput1,iInput2},'lags','ccf','ccfBounds',...
%         'bootstrapCcf','bootstrapCcfBounds','-append');
%     end
% end

disp('Finished calculating partial correlation...')

