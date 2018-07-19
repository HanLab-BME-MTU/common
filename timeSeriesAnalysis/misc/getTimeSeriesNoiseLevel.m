function  outTS = getTimeSeriesNoiseLevel(TS,varargin)

%This function estimates the local noise level (SNR) in a Time Series
%
%WARNING: Noise here is represented by the fastest mode of the signal.
%         Significant changes in amplitude is taken as assumption to real identify signal at the fastest frequency              
%
%Usage:
%       [outTS] = getTimeSeriesNoiseLevel(TS,varargin)
%
%Input:
%       TS        - vector with time observations
%       alpha     - alpha level for testing the IMF's (see testImf)
%       winSize   - size of the sliding window used to average the local SNR
%       noiseStd  - number of std of the local noise distribution that will be used to define the local noise limits(Default=1)
%       plotYes   - if 1, plots the time series with the noise level
%       frameRate - just to set the right time axis in the plot
%
%Output:
%       outTS.Up          - Upper noise limit (same size as TS)
%       outTS.Dw          - Lower noise limit (same size as TS)
%       outTS.localSNR    - local SNR(same size as TS)
%       outTS.globalSNR   - average SNR for the entire tsMap (scalar)
%
% See also: getSNRmap, testImf
%Marco Vilela, 2012



ip=inputParser;
ip.addRequired('TS',@(x) isnumeric(x));
ip.addParamValue('alpha',   .05,@isscalar);
ip.addParamValue('winSize', 10,@isscalar);
ip.addParamValue('noiseStd', 1,@isscalar);
ip.addParamValue('plotYes', 0,@isscalar);
ip.addParamValue('frameRate',1,@isscalar);

ip.parse(TS,varargin{:});
alpha     = ip.Results.alpha;
winSize   = ip.Results.winSize;
noiseStd  = ip.Results.noiseStd;
plotYes   = ip.Results.plotYes;
frameRate = ip.Results.frameRate;

if mod(winSize,2) == 0
    
    winSize = winSize + 1;%Odd number
    
end



TS       = TS(:)';
nPoint   = numel(TS);
outTS.Up = nan(size(TS));
outTS.Dw = nan(size(TS));


if ~isempty(TS)
    
    imf       = emd(TS);
    Mu        = mean(imf(1,:));%Fastest mode
    [~,noise] = testImf(imf,'alpha',alpha);
    
    %Local variance calculation
    
    localSignalVar  = slidingWindowFilter(TS',winSize,@(x) nanvar(x,1,2));
    localNoiseVar   = slidingWindowFilter(noise,winSize,@(x) nanvar(x,1,2));
    outTS.localSNR  = localSignalVar./localNoiseVar;
    outTS.globalSNR = nanvar(TS)/nanvar(noise);
    
    %Defining the lower and upper noise confidence bands
    outTS.Up  = Mu + sqrt(localNoiseVar)*noiseStd;
    outTS.Dw  = Mu - sqrt(localNoiseVar)*noiseStd;
    %*************************
    if plotYes
                
        alphaT    = 1;
        figure
        xAxis = frameRate*( 1:nPoint );
        plot(xAxis,TS,'b','LineWidth',2)
        xlabel('Time [sec]')
        ylabel('Signal [A.U]')
        hold on
        
        
        fill([xAxis fliplr(xAxis)]',[outTS.Up;flipud(outTS.Dw)],'r','FaceAlpha', alphaT,'linestyle','none');
        disp(['If you are going to import this figure into Illustrator, change the alphaT to 1' ...
              'DO NOT CHECK IT IN'])
        %The shading can be done once the figure is imported into illustrator by altering the transparency of the red area
        
        h = legend({'Signal','Noise Level'});
        legend(h,'boxoff')
        
    end
    
end
