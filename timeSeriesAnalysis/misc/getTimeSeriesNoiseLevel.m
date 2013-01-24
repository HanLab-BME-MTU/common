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
%       outLevel  - Std for the outliear detection    (see detectOutliers)
%       trend     - type of trend to be removed       (see removeMean) 
%       minLen    - minimal time series length
%       winSize   - size of the sliding window used to average the local SNR
%       noiseStd  - number of std of the local noise distribution that will be used to define the local noise limits
%       plotYes   - if 1, plots the time series with the noise level
%       frameRate - just to set the right time axis in the plot
%
%Output:
%       outTS.points   - processed TS (trend, mean, nan removal) - same size as TS
%       outTS.interval - 
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
ip.addParamValue('outLevel',7,@isscalar);
ip.addParamValue('trendType',   1,@isscalar);
ip.addParamValue('minLength',  30,@isscalar);
ip.addParamValue('winSize', 10,@isscalar);
ip.addParamValue('noiseStd', 1,@isscalar);
ip.addParamValue('plotYes', 0,@isscalar);
ip.addParamValue('frameRate',1,@isscalar);

ip.parse(TS,varargin{:});
alpha     = ip.Results.alpha;
outLevel  = ip.Results.outLevel;
trend     = ip.Results.trendType;
minLen    = ip.Results.minLength;
winSize   = ip.Results.winSize;
noiseStd  = ip.Results.noiseStd;
plotYes   = ip.Results.plotYes;
frameRate = ip.Results.frameRate;

if mod(winSize,2) == 0
    
    winSize = winSize + 1;%Odd number
    
end


    % IT SHOULD NOT HAVE PRE-PROCESSING IN IT

TS             = TS(:)';
outTS.points   = timeSeriesPreProcessing(TS,'gapSize',1,'trendType',trend);
outTS.interval = 1:length(TS);
%TS(detectOutliers(TS,outLevel)) = NaN;
%[outTS.points,outTS.interval]   = removeMeanTrendNaN(TS,'trendType',trend,'minLength',minLen);
%outTS.points                    = cell2mat(outTS.points);
%outTS.interval                  = cell2mat(outTS.interval);

outTS.Up                        = nan(size(TS));
outTS.Dw                        = nan(size(TS));
TSwork                          = outTS.points';

if ~isempty(TSwork)
    
    imf       = emd(TSwork);
    Mu        = mean(imf(1,:));
    [~,noise] = testImf(imf,'alpha',alpha);
    
    %Local variance calculation
    
    localSignalVar  = slidingWindowFilter(TSwork,winSize,@(x) nanvar(x,1,2));
    localNoiseVar   = slidingWindowFilter(noise,winSize,@(x) nanvar(x,1,2));
    outTS.localSNR  = localSignalVar./localNoiseVar;
    outTS.globalSNR = nanvar(TSwork)/nanvar(noise);
    
    %Defining the lower and upper noise confidence bands
    outTS.Up(outTS.interval)  = Mu + sqrt(localNoiseVar)*noiseStd;
    outTS.Dw(outTS.interval)  = Mu - sqrt(localNoiseVar)*noiseStd;
    %*************************
    if plotYes
                
        alphaT    = 0.2;
        figure
        xAxis = frameRate*( outTS.interval );
        plot(xAxis,TSwork,'LineWidth',2)
        xlabel('Time [sec]')
        ylabel('Signal [A.U]')
        hold on
        
        
        fill([xAxis fliplr(xAxis)],[outTS.Up flipud(outTS.Dw)],'r','FaceAlpha', alphaT,'linestyle','none');
        disp(['If you are going to import this figure into Illustrator, comment the above and uncomment the line below' ...
              'DO NOT FORGET TO REVERSE THE COMMENTS BEFORE svn ci'])
        %The shading can be done once the figure is imported into illustrator by altering the transparency of the red area
        
        %fill([xAxis;flipud(xAxis)],[Up;flipud(Dw)],'r','linestyle','none');
        h = legend({'Detrended Signal','Noise Level'});
        legend(h,'boxoff')
        
    end
    
end
