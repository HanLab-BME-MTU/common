function [meanAll]=errorBarPlotCellArray(cellArrayData,nameList,convertFactor,varargin)
% function [meanAll]=errorBarPlotCellArray(FAarea) automatically converts cell array
% format input to matrix input to use matlab function 'errorbar'
% input: cellArrayData      cell array data
%           nameList        cell array containing name of each
%                             condition (eg: {'condition1' 'condition2' 'condition3'})
%                           This can be also numeric data, in case of which
%                           each value will be used for x-axis.
%           convertFactor     conversion factor for physical unit (ex.
%           pixelSize, timeInterval etc...)
% Sangyoon Han, March 2016
ip = inputParser;
ip.CaseSensitive = false;
addRequired(ip,'cellArrayData');
addOptional(ip,'nameList',arrayfun(@(x) num2str(x),(1:numel(cellArrayData)),'UniformOutput',false));
addOptional(ip,'convertFactor',1);
addParameter(ip,'colorful',true);
addParameter(ip,'color','r');
addParameter(ip,'Stat',true);
addParameter(ip,'markerSize',6);
addParameter(ip,'ax',gca);
addParameter(ip,'horizontalPlot',false);
addParameter(ip,'usePairedTTest',false);

parse(ip,cellArrayData,nameList,convertFactor,varargin{:});
MSize = ip.Results.markerSize;
colorful = ip.Results.colorful;
colorSelected = ip.Results.color;
doStat = ip.Results.Stat;
usePairedTTest = ip.Results.usePairedTTest;

[lengthLongest]=max(cellfun(@(x) length(x),cellArrayData));
numConditions = numel(cellArrayData);
matrixData = NaN(lengthLongest,numConditions);
for k=1:numConditions
    matrixData(1:length(cellArrayData{k}),k) = cellArrayData{k};
end

if numConditions<18
    colorSwitch = 'qkjlrabnfilkpuvdt';
    colors = extendedColors(colorSwitch);
else
    colors = distinguishable_colors(numConditions,'w');
end

matrixData=matrixData*convertFactor;
meanAll = nanmean(matrixData);

if colorful
    for ii=1:numConditions
        if iscell(nameList)
            errorbar(ii,nanmean(matrixData(:,ii)),...
                stdErrMean(matrixData(:,ii)),...
                'Marker', 'x','MarkerSize',MSize,'Color', colors(ii,:));
        else
            errorbar(nameList(ii),nanmean(matrixData(:,ii)),...
                stdErrMean(matrixData(:,ii)),...
                'Marker', 'x','MarkerSize',MSize,'Color', colors(ii,:));
        end
        hold on
    end
else
    if iscell(nameList)
        errorbar(1:numConditions,nanmean(matrixData),...
            stdErrMean(matrixData),...
            'Marker', 'x','MarkerSize',MSize,'Color', colorSelected);
    else
        errorbar(nameList,nanmean(matrixData),...
            stdErrMean(matrixData),...
            'Marker', 'x','MarkerSize',MSize,'Color', colorSelected);
    end
end

set(gca,'XTickLabelRotation',45)
if iscell(nameList)
    set(gca,'XTickLabel',nameList)
    set(gca,'XTick',1:numel(cellArrayData))
    xlim([0.5 numConditions+0.5])
else
    set(gca,'XTickLabel',nameList)
    set(gca,'XTick',nameList)
    xlim([min(0,nameList(1)) nameList(end)*1.1])
end

hold on
% perform ranksum test for every single combination
maxPoint =cellfun(@nanmean,cellArrayData)+cellfun(@(x) stdErrMean(x),cellArrayData); %/sqrt(length(x))
maxPoint2 = nanmax(maxPoint(:));
lineGap=maxPoint2*0.03;
qUsed=[];
xUsed=[];
xGap = 0.02;
ax = gca;
if doStat
    if iscell(nameList)
        xi = 1:numel(nameList);
    else
        xi = nameList;
    end
    for k=1:(numConditions-1)
        q=-2*lineGap;
        for ii=k+1:numConditions
            if numel(cellArrayData{k})>1 && numel(cellArrayData{ii})>1
                if usePairedTTest
                    method = 'paired t-test';
                    [~,p] = ttest(cellArrayData{k},cellArrayData{ii});
                elseif kstest(cellArrayData{k}) % this means the test rejects the null hypothesis
                    method = 'ranksum test';
                    [p]=ranksum(cellArrayData{k},cellArrayData{ii});
                else
                    method = 'unpaired t-test';
                    [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
                end
                if (p<0.05) 
                    q = max(cellfun(@(x) nanmean(x)+stdErrMean(x),cellArrayData(k:ii)));
                    tol=0.95*lineGap;

                    qOverlap=abs(qUsed-q)<tol;
                    xOverlap=false(1,size(xUsed,1));
                    for kk=1:size(xUsed,1)
                        xOverlap(kk)=(k<xUsed(kk,2) & ii>xUsed(kk,1));
                    end
                    a = find(qOverlap & xOverlap,1);
                    while ~isempty(a) %ismember(q,qUsed) 
                        q=qUsed(a)+lineGap;
                        qOverlap=abs(qUsed-q)<tol;
                        a = find(qOverlap & xOverlap,1);
                    end
                    qUsed = [qUsed q];
                    xUsed = [xUsed; k ii];
                    q=q+lineGap*0.1;
                    line(ax,[k+xGap ii-xGap], ones(1,2)*(q),'Color','k')    
                    q=q+lineGap*0.3;
%                         text(ax,floor((k+ii)/2)+0.3, maxPoint2+q,['p=' num2str(p,'%2.2e') ' (r)'])
                    text(ax,((k+ii)/2)-0.4, q,['p=' num2str(p,2)]) %'%2.2e'
%                         line(ax,[xi(k)+xGap xi(ii)-xGap], ones(1,2)*(maxPoint2+q),'Color','k')    
%                         q=q+lineGap;
%                         text(ax,floor((xi(k)+xi(ii))/2)+0.3, maxPoint2+q,['p=' num2str(p,'%2.2e') ' (r)'])
                end

            end
        end
    end
else
    q=0;
end
set(gca,'FontSize',6)
set(findobj(gca,'Type','Text'),'FontSize',5)

q=q+lineGap*3;
minPoint =cellfun(@nanmedian,cellArrayData)-cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
minPoint2 = nanmin(minPoint(:));

yMin = max(minPoint2-lineGap*2, 0);
ylim([yMin maxPoint2+q+lineGap*8])
text(ax,ii-1,yMin+lineGap,method)

