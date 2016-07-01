function []=boxPlotCellArray(cellArrayData,nameList,convertFactor,notchOn)
% function []=boxPlotCellArray(FAarea) automatically converts cell array
% format input to matrix input to use matlab function 'boxplot'
% input: cellArrayData      cell array data
%           nameList            cell array containing name of each
%                                       condition (ex: {'condition1' 'condition2' 'condition3'})
%           convertFactor     conversion factor for physical unit (ex.
%           pixelSize, timeInterval etc...)
% Sangyoon Han, March 2016
[lengthLongest]=max(cellfun(@(x) length(x),cellArrayData));
numConditions = numel(cellArrayData);
matrixData = NaN(lengthLongest,numConditions);
for k=1:numConditions
    matrixData(1:length(cellArrayData{k}),k) = cellArrayData{k};
end
if nargin<4
    notchOn=true;
end
if nargin<3
    convertFactor = 1;
    notchOn=true;
end
if nargin<2
    nameList=arrayfun(@(x) num2str(x),(1:numConditions),'UniformOutput',false);
    convertFactor = 1;
    notchOn=true;
end
boxWidth=0.5;
whiskerRatio=1;
matrixData=matrixData*convertFactor;
if notchOn %min(sum(~isnan(matrixData),1))>20 || 
    boxplot(matrixData,'whisker',whiskerRatio,'notch','on',...
        'labels',nameList,'symbol','','widths',boxWidth,'jitter',1,'colors','k', 'labelorientation','inline')
else % if the data is too small, don't use notch
    boxplot(matrixData,'whisker',whiskerRatio*0.5,'notch','off',...
        'labels',nameList,'symbol','','widths',boxWidth,'jitter',1,'colors','k', 'labelorientation','inline')
end
    
set(findobj(gca,'LineStyle','--'),'LineStyle','-')
set(findobj(gca,'tag','Median'),'LineWidth',2)

hold on
% perform ranksum test for every single combination
maxPoint = quantile(matrixData,[0.25 0.75]);
maxPoint2 = maxPoint(2,:)+(maxPoint(2,:)-maxPoint(1,:))*whiskerRatio;
maxPoint2 = max(maxPoint2);
lineGap=maxPoint2*0.05;
q=0;
for k=1:(numConditions-1)
    for ii=k+1:numConditions
        q=q+lineGap;
        line([k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
        q=q+lineGap;
        if kstest(cellArrayData{k}) % this means the test rejects the null hypothesis
            [p]=ranksum(cellArrayData{k},cellArrayData{ii});
            text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) ' (r)'])
       else
            [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
            text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) ' (t)'])
       end
    end
end
q=q+lineGap*3;
ylim([0 maxPoint2+q])

