function []=barPlotCellArray(cellArrayData,nameList,convertFactor)
% function []=barPlotCellArray(cellArrayData,nameList,convertFactor) automatically converts cell array
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
errorbar(1:numConditions,cellfun(@nanmean,cellArrayData)*convertFactor,cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData)*convertFactor,'Marker','none','LineStyle','none','Color','k','LineWidth',1.5); %SEM option 
% errorbar(1:numConditions, cellfun(@nanmean,cellArrayData)*convertFactor,cellfun(@(x) nanstd(x),cellArrayData)*convertFactor,'Marker','none','LineStyle','none','Color','k','LineWidth',1.5); %SD option 
try
    nameListNew = cellfun(@(x,y) [x '(N=' num2str(length(y)) ')'],nameList,cellArrayData,'UniformOutput', false);
catch
    nameListNew = cellfun(@(x,y) [x '(N=' num2str(length(y)) ')'],nameList,cellArrayData','UniformOutput', false);
end
set(gca,'XTickLabel',nameListNew)
set(gca,'XTickLabelRotation',45)

%% perform ranksum test for every single combination
meanPoints =cellfun(@nanmean,cellArrayData);
maxPoint = max(meanPoints);
minPoint = min(meanPoints);
if maxPoint>0 && minPoint>0
    maxPointWithStd =cellfun(@nanmean,cellArrayData)+cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
    extremePoint2 = nanmax(maxPointWithStd(:))*convertFactor;
elseif maxPoint<0 && minPoint<0
    maxPointWithStd =cellfun(@nanmean,cellArrayData)-cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
    extremePoint2 = nanmin(maxPointWithStd(:))*convertFactor;
else
    if maxPoint>abs(minPoint)
        maxPointWithStd =cellfun(@nanmean,cellArrayData)+cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
        extremePoint2 = nanmax(maxPointWithStd(:))*convertFactor;
    else
        maxPointWithStd =cellfun(@nanmean,cellArrayData)-cellfun(@(x) nanstd(x)/sqrt(length(x)),cellArrayData);
        extremePoint2 = nanmin(maxPointWithStd(:))*convertFactor;
    end
end
lineGap=abs(extremePoint2*0.05);
q=0;
for k=1:(numConditions-1)
    for ii=k+1:numConditions
        if numel(cellArrayData{k})>1 && numel(cellArrayData{ii})>1
            if kstest(cellArrayData{k}) % this means the test rejects the null hypothesis
                [p]=ranksum(cellArrayData{k},cellArrayData{ii});
                if p<0.05
                    if extremePoint2>0
                        q=q+lineGap;
                        line([k ii], ones(1,2)*(extremePoint2+q),'Color','k')    
                        q=q+lineGap;
                        text(floor((k+ii)/2), extremePoint2+q,['p=' num2str(p) ' (r)'])
                    else
                        q=q+lineGap;
                        line([k ii], ones(1,2)*(extremePoint2-q),'Color','k')    
                        q=q+lineGap;
                        text(floor((k+ii)/2), extremePoint2-q,['p=' num2str(p) ' (r)'])
                    end
                end
            else
                [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
                if p<0.05
                    if extremePoint2>0
                        q=q+lineGap;
                        line([k ii], ones(1,2)*(extremePoint2+q),'Color','k')    
                        q=q+lineGap;
                        text(floor((k+ii)/2), extremePoint2+q,['p=' num2str(p) ' (t)'])
                    else
                        q=q+lineGap;
                        line([k ii], ones(1,2)*(extremePoint2-q),'Color','k')    
                        q=q+lineGap;
                        text(floor((k+ii)/2), extremePoint2-q,['p=' num2str(p) ' (t)'])
                    end
                end
            end
        end
    end
end
% for k=1:(numConditions-1)
%     for ii=k+1:numConditions
%         q=q+lineGap;
%         line([k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
%         q=q+lineGap;
%         if ~kstest(cellArrayData{k})
%             [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
%             text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) '(t)'])
%        else
%             [p]=ranksum(cellArrayData{k},cellArrayData{ii});
%             text(floor((k+ii)/2), maxPoint2+q,['p=' num2str(p) '(r)'])
%        end
%     end
% end
q=q+lineGap*3;
if extremePoint2>0
    if minPoint>0
        ylim([0 extremePoint2+q])
    else
        ylim([minPoint-lineGap extremePoint2+q])
    end
else
    if minPoint<0
        ylim([extremePoint2-q 0])
    else
        ylim([extremePoint2-q minPoint+lineGap])
    end
end

