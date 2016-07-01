function []=barPlotCellArray(cellArrayData,nameList,convertFactor)
% function []=barPlotCellArray(FAarea) automatically converts cell array
% format input to matrix input to use matlab function 'boxplot'
% input: cellArrayData      cell array data
%           nameList            cell array containing name of each
%                                       condition (ex: {'condition1' 'condition2' 'condition3'})
%           convertFactor     conversion factor for physical unit (ex.
%           pixelSize, timeInterval etc...)
% Sangyoon Han, March 2016
if nargin<3
    convertFactor = 1;
end
numConditions=numel(cellArrayData);
if nargin<2
    nameList=arrayfun(@(x) num2str(x),(1:numConditions),'UniformOutput',false);
end
bar(1:numConditions, cellfun(@nanmean,cellArrayData)*convertFactor,'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',1.5); hold on
errorbar(1:numConditions, cellfun(@nanmean,cellArrayData)*convertFactor,cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData)*convertFactor,'Marker','none','LineStyle','none','Color','k','LineWidth',1.5);
set(gca,'XTickLabel',nameList)

%% perform ranksum test for every single combination
maxPoint =cellfun(@nanmean,cellArrayData)+cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
maxPoint2 = nanmax(maxPoint(:));
lineGap=maxPoint2*0.05;
q=0;
for k=1:(numConditions-1)
    for ii=k+1:numConditions
        q=q+lineGap;
        line([k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
        q=q+lineGap;
        if ~kstest(cellArrayData{k})
            [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
            text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) '(t)'])
       else
            [p]=ranksum(cellArrayData{k},cellArrayData{ii});
            text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) '(r)'])
       end
    end
end
q=q+lineGap*3;
ylim([0 maxPoint2+q])

