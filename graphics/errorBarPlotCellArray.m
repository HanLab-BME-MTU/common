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
addParameter(ip,'markerSize',6);
addParameter(ip,'ax',gca);
addParameter(ip,'horizontalPlot',false);
parse(ip,cellArrayData,nameList,convertFactor,varargin{:});
MSize = ip.Results.markerSize;

[lengthLongest]=max(cellfun(@(x) length(x),cellArrayData));
numConditions = numel(cellArrayData);
matrixData = NaN(lengthLongest,numConditions);
for k=1:numConditions
    matrixData(1:length(cellArrayData{k}),k) = cellArrayData{k};
end

colorSwitch = 'qkjlrabnfilkpuvdt';
color = extendedColors(colorSwitch);

matrixData=matrixData*convertFactor;
meanAll = nanmean(matrixData);

for ii=1:numConditions
    if iscell(nameList)
        errorbar(ii,nanmean(matrixData(:,ii)),...
            stdErrMean(matrixData(:,ii)),...
            'Marker', 'x','MarkerSize',MSize,'Color', color(ii,:));
    else
        errorbar(nameList(ii),nanmean(matrixData(:,ii)),...
            stdErrMean(matrixData(:,ii)),...
            'Marker', 'x','MarkerSize',MSize,'Color', color(ii,:));
    end
    hold on
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
maxPoint =cellfun(@nanmedian,cellArrayData)+cellfun(@(x) nanstd(x),cellArrayData); %/sqrt(length(x))
maxPoint2 = nanmax(maxPoint(:));
lineGap=maxPoint2*0.03;
xGap = 0.02;
ax = gca;
if iscell(nameList)
    xi = 1:numel(nameList);
else
    xi = nameList;
end
for k=1:(numConditions-1)
    q=-2*lineGap;
    for ii=k+1:numConditions
        if numel(cellArrayData{k})>1 && numel(cellArrayData{ii})>1
            if kstest(cellArrayData{k}) % this means the test rejects the null hypothesis
                [p]=ranksum(cellArrayData{k},cellArrayData{ii});
                if (p<0.05) 
                    q=q+lineGap;
                    line(ax,[xi(k)+xGap xi(ii)-xGap], ones(1,2)*(maxPoint2+q),'Color','k')    
                    q=q+lineGap;
                    text(ax,floor((xi(k)+xi(ii))/2)+0.3, maxPoint2+q,['p=' num2str(p,'%2.2e') ' (r)'])
                end
            else
                [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
                if (p<0.05 ) 
                    q=q+lineGap;
                    line(ax,[xi(k)+xGap xi(ii)-xGap], ones(1,2)*(maxPoint2+q),'Color','k')    
                    q=q+lineGap;
                    text(ax,floor((xi(k)+xi(ii))/2)+0.3, maxPoint2+q,['p=' num2str(p,'%3.2e') ' (t)'])
                end
            end
        end
    end
end

set(gca,'FontSize',6)
set(findobj(gca,'Type','Text'),'FontSize',5)

q=q+lineGap*3;
minPoint =cellfun(@nanmedian,cellArrayData)-cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
minPoint2 = nanmin(minPoint(:));

yMin = min(minPoint2-lineGap*2, 0);
ylim([yMin maxPoint2+q+lineGap*8])
