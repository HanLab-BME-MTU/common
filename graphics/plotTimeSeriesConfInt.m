function plotTimeSeriesConfInt(xSeries,curArray, varargin)
% function plotTimeSeriesConfInt(xSeries,curArray, varargin) plots mean
% time series of curArray (nanmean(curArray,1)) with confidence interval
% with a shade.
% input
%   xSeries:    1xn vector for time
%   curArray:   mxn matrix containing m different time sereis with n
%               maximum 
% example:
%   plotTimeSeriesConfInt((1:nSampleFrames)*tInterval,forceArray, 
%                         'Color', [240/255 128/255 128/255])

ip = inputParser;
ip.addParameter('Color',[0.5 0.5 0.5],@isnumeric);
% ip.addParamValue('YLim',[],@isnumeric || @isempty);
% ip.addParamValue('tInterval',1,@isnumeric);

ip.parse(varargin{:});

colorUsed = ip.Results.Color;
% YLim = ip.Results.YLim;
% tInterval = ip.Results.tInterval;

ciColor = colorUsed*1.5;%[153/255 255/255 51/255];
% Saturate
for ii=1:3
    if ciColor(ii)>1
        ciColor(ii)=1;
    end
end
meanColor = colorUsed*0.5;%[0/255 102/255 0];
curLT =length(xSeries);
tInterval = xSeries(2)-xSeries(1);
normalizedXseries=xSeries/tInterval;

curMeanSig = nanmean(curArray,1);
curSEM = nanstd(curArray,1)/sqrt(size(curArray,1));
curTScore = tinv([0.025 0.975],size(curArray,1)-1);
curCI_upper = curMeanSig + curTScore*curSEM;
curCI_lower = curMeanSig - curTScore*curSEM;
fill([normalizedXseries(1):normalizedXseries(end) normalizedXseries(end):-1:normalizedXseries(1)]*tInterval,[curCI_upper fliplr(curCI_lower)],ciColor,'EdgeColor',ciColor),hold on
plot((normalizedXseries(1):normalizedXseries(end))*tInterval,curMeanSig,'Linewidth',2,'Color',meanColor)





