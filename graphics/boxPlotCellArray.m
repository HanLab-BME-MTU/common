function curAxes=boxPlotCellArray(cellArrayData,nameList,convertFactor,notchOn,plotIndivPoint,forceShowP)
% function []=boxPlotCellArray(FAarea) automatically converts cell array
% format input to matrix input to use matlab function 'boxplot'
% input: cellArrayData      cell array data
%           nameList            cell array containing name of each
%                                       condition (ex: {'condition1' 'condition2' 'condition3'})
%           convertFactor     conversion factor for physical unit (ex.
%           pixelSize, timeInterval etc...)
% Sangyoon Han, March 2016
[lengthLongest]=max(cellfun(@(x) length(x),cellArrayData));
%If there is no data, exclude them in the plot
idEmptyData = cellfun(@isempty,cellArrayData);
cellArrayData(idEmptyData)=[];
nameList(idEmptyData)=[];
numConditions = numel(cellArrayData);
matrixData = NaN(lengthLongest,numConditions);
for k=1:numConditions
    matrixData(1:length(cellArrayData{k}),k) = cellArrayData{k};
end
if isempty(matrixData)
    text(0.4,0.5,'Empty Data')
    return
end
if nargin<6
    forceShowP=false;
end
if nargin<5
    forceShowP=false;
    plotIndivPoint = true;
end
if nargin<4
    forceShowP=false;
    notchOn=true;
    plotIndivPoint = true;
end
if nargin<3
    forceShowP=false;
    convertFactor = 1;
    notchOn=true;
end
if nargin<2
    nameList=arrayfun(@(x) num2str(x),(1:numConditions),'UniformOutput',false);
    convertFactor = 1;
    notchOn=true;
    forceShowP=false;
end
boxWidth=0.5;
whiskerRatio=1.5;
matrixData=matrixData*convertFactor;
nameListNew = cellfun(@(x,y) [x '(N=' num2str(length(y)) ')'],nameList,cellArrayData,'UniformOutput', false);

% individual data jitter plot
numCategories = size(matrixData,2);
% width=boxWidth/2;
for ii=1:numCategories
    xData = ii+0.1*boxWidth*(randn(size(matrixData,1),1));
    scatter(xData,matrixData(:,ii),'filled','MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','SizeData',2)
    hold on
end
alpha(.5)

if notchOn %min(sum(~isnan(matrixData),1))>20 || 
    boxplot(matrixData,'whisker',whiskerRatio,'notch','on',...
        'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',1,'colors','k');%, 'labelorientation','inline');
else % if the data is too small, don't use notch
    boxplot(matrixData,'whisker',whiskerRatio*0.5,'notch','off',...
        'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',0,'colors','k');%, 'labelorientation','inline');
end
    
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)
% set(gca,'XTick',1:numel(nameList))
% set(gca,'XTickLabel',nameList)
set(gca,'XTickLabelRotation',45)
curAxes = gca;

% hold on
% perform ranksum test for every single combination
maxPoint = quantile(matrixData,[0.25 0.75]);
if size(matrixData,2)>1
    maxPoint2 = maxPoint(2,:)+(maxPoint(2,:)-maxPoint(1,:))*whiskerRatio;
    maxPoint2 = max(maxPoint2);
else
    maxPoint2 = maxPoint(2)+(maxPoint(2)-maxPoint(1))*whiskerRatio;
end
lineGap=maxPoint2*0.05;
q=0;
for k=1:(numConditions-1)
    for ii=k+1:numConditions
        if numel(cellArrayData{k})>1 && numel(cellArrayData{ii})>1
            if kstest(cellArrayData{k}) % this means the test rejects the null hypothesis
                [p]=ranksum(cellArrayData{k},cellArrayData{ii});
                if p<0.05 || forceShowP
                    q=q+lineGap;
                    line([k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
                    q=q+lineGap;
                    text(floor((k+ii)/2)+0.3, maxPoint2+q,['p=' num2str(p) ' (r)'])
                end
            else
                [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
                if p<0.05 || forceShowP
                    q=q+lineGap;
                    line([k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
                    q=q+lineGap;
                    text(floor((k+ii)/2)+0.3, maxPoint2+q,['p=' num2str(p) ' (t)'])
                end
            end
        end
    end
end
q=q+lineGap*3;
minPoint = quantile(matrixData,[0.25 0.75]);
if size(matrixData,2)>1
    minPoint2 = minPoint(1,:)-(minPoint(2,:)-minPoint(1,:))*whiskerRatio;
    minPoint2 = min(minPoint2);
else
    minPoint2 = minPoint(2)-(minPoint(2)-minPoint(1))*whiskerRatio;
end
if plotIndivPoint && ~forceShowP
    maxPoint2 = quantile(matrixData(:),0.99);
else
    maxPoint2 = maxPoint2+q;
end
if maxPoint2>minPoint2-lineGap*2
    ylim([minPoint2-lineGap*2 maxPoint2])
else
    ylim auto
end

set(gca,'FontSize',7)
set(findobj(gca,'Type','Text'),'FontSize',6)

