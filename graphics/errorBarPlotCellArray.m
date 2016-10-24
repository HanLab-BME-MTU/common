function []=errorBarPlotCellArray(cellArrayData,nameList,convertFactor)
% function []=errorBarPlotCellArray(FAarea) automatically converts cell array
% format input to matrix input to use matlab function 'errorbar'
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
if nargin<3
    convertFactor = 1;
end
if nargin<2
    nameList=arrayfun(@(x) num2str(x),(1:numConditions),'UniformOutput',false);
    convertFactor = 1;
end
matrixData=matrixData*convertFactor;
errorbar(nanmedian(matrixData),...
    stdErrMean(matrixData),'rx')
set(gca,'XTickLabel',nameList)
set(gca,'XTick',1:numel(cellArrayData))
set(gca,'XTickLabelRotation',45)

hold on
% perform ranksum test for every single combination
maxPoint =cellfun(@nanmedian,cellArrayData)+cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
maxPoint2 = nanmax(maxPoint(:));
lineGap=maxPoint2*0.05;
q=0;
for k=1:(numConditions-1)
    for ii=k+1:numConditions
        if numel(cellArrayData{k})>1 && numel(cellArrayData{ii})>1
            if kstest(cellArrayData{k}) % this means the test rejects the null hypothesis
                [p]=ranksum(cellArrayData{k},cellArrayData{ii});
                if p<0.05
                    q=q+lineGap;
                    line([k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
                    q=q+lineGap;
                    text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) ' (r)'])
                end
            else
                [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
                if p<0.05
                    q=q+lineGap;
                    line([k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
                    q=q+lineGap;
                    text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) ' (t)'])
                end
            end
        end
    end
end
set(gca,'FontSize',6)
set(findobj(gca,'Type','Text'),'FontSize',4)

q=q+lineGap*3;
minPoint =cellfun(@nanmedian,cellArrayData)-cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
minPoint2 = nanmin(minPoint(:));

ylim([minPoint2-lineGap*2 maxPoint2+q])

