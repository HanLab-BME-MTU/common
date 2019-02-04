function [sucess]=boxPlotCellArray(cellArrayData,nameList,convertFactor,notchOn,plotIndivPoint,forceShowP,markerSize,varargin)
% function []=boxPlotCellArray(cellArrayData,nameList,convertFactor,notchOn,plotIndivPoint,forceShowP) automatically converts cell array
% format input to matrix input to use matlab function 'boxplot'
% input: cellArrayData          cell array data
%           nameList            cell array containing name of each
%                               condition (ex: {'condition1' 'condition2' 'condition3'})
%           convertFactor       conversion factor for physical unit (ex.
%                               pixelSize, timeInterval etc...)
%           forceShowP          0 if you want to show only significant p
%                               1 if you want to show all p 
%                               2 if you do not want to show any p
% Sangyoon Han, March 2016
[lengthLongest]=max(cellfun(@(x) length(x),cellArrayData));
%If there is no data, exclude them in the plot
idEmptyData = cellfun(@isempty,cellArrayData);
cellArrayData(idEmptyData)=[];
sucess=false;
% Return when one of columns are all NaNs
isAnyGroupAllNaNs = cellfun(@(x) all(isnan(x)), cellArrayData);
if any(isAnyGroupAllNaNs)
    % Making the 1xN vector if isAnyGroupAllNaNs is Nx1 vector
    if size(isAnyGroupAllNaNs,1)>1
        isAnyGroupAllNaNs=isAnyGroupAllNaNs';
    end
    disp([num2str(find(isAnyGroupAllNaNs)) 'th group contains all NaNs. Returning...'])
    return
end

numConditions = numel(cellArrayData);
matrixData = NaN(lengthLongest,numConditions);
for k=1:numConditions
    matrixData(1:length(cellArrayData{k}),k) = cellArrayData{k};
end
if all(isnan(matrixData(:)))
    disp('All data are NaNs. Returning...')
    return
end

ip = inputParser;
ip.CaseSensitive = false;
addRequired(ip,'cellArrayData');
addOptional(ip,'nameList',arrayfun(@(x) num2str(x),(1:numConditions),'UniformOutput',false));
addOptional(ip,'convertFactor',1);
addOptional(ip,'notchOn',true);
addOptional(ip,'plotIndivPoint',true);
addOptional(ip,'forceShowP',false);
addOptional(ip,'markerSize',2);
addParameter(ip,'ax',gca);
addParameter(ip,'horizontalPlot',false);
parse(ip,cellArrayData,nameList,convertFactor,notchOn,plotIndivPoint,forceShowP,markerSize,varargin{:});
ax=ip.Results.ax;
horizontalPlot=ip.Results.horizontalPlot;
markerSize = ip.Results.markerSize;
forceShowP = ip.Results.forceShowP;
notchOn = ip.Results.notchOn;
plotIndivPoint = ip.Results.plotIndivPoint;
convertFactor = ip.Results.convertFactor;
% if nargin<7
%     markerSize=2;
% end
% if isempty(matrixData)
%     text(0.4,0.5,'Empty Data')
%     return
% end
% if nargin<6
%     markerSize=2;
%     forceShowP=false;
% end
% if nargin<5
%     markerSize=2;
%     forceShowP=false;
%     plotIndivPoint = true;
% end
% if nargin<4
%     markerSize=2;
%     forceShowP=false;
%     notchOn=true;
%     plotIndivPoint = true;
% end
% if nargin<3
%     markerSize=2;
%     forceShowP=false;
%     convertFactor = 1;
%     notchOn=true;
% end
% if nargin<2
%     markerSize=2;
%     nameList=arrayfun(@(x) num2str(x),(1:numConditions),'UniformOutput',false);
%     convertFactor = 1;
%     notchOn=true;
%     forceShowP=false;
% end
nameList(idEmptyData)=[];

boxWidth=0.5;
scatterWidth=boxWidth*2;
whiskerRatio=1.5;
matrixData=matrixData*convertFactor;
nameListNew = cellfun(@(x,y) [x '(N=' num2str(sum(~isnan(y))) ')'],nameList,cellArrayData,'UniformOutput', false);

onlyOneDataAllGroups=false;
numCategories = numel(cellArrayData);
onlyOneDataPerEachGroup=false(1,numCategories);
if plotIndivPoint
    % individual data jitter plot
%     numCategories = size(matrixData,2);
    % width=boxWidth/2;
    for ii=1:numCategories
        Nall{ii} = histcounts(matrixData(:,ii));
        uArrayAll{ii} = unique(matrixData(~isnan(matrixData(:,ii)),ii));
    end
    if mean(cellfun(@(x,y) length(y)/length(x),uArrayAll,cellArrayData))>3
        Nmax = max(cellfun(@max,cellfun(@(x,y) arrayfun(@(z) sum(y==z),x) ,uArrayAll,cellArrayData,'unif',false)));
    else
        Nmax = max(cellfun(@max,Nall));
    end
    for ii=1:numCategories
        curNumData = numel(matrixData(:,ii));
        if curNumData==1
            plot(ii,matrixData(:,ii),'k.')
            onlyOneDataPerEachGroup(ii)=true;
        else
            xData = ii+0.1*boxWidth*(randn(size(matrixData,1),1));
            % Need to take care of xData more wisely
            % Going with matrixData(:,ii)
            uCurArray = unique(matrixData(~isnan(matrixData(:,ii)),ii));
            if 3*length(uCurArray) < sum(~isnan(matrixData(:,ii)))
                N = arrayfun(@(x) sum(matrixData(:,ii)==x),uCurArray)';
                edges=[uCurArray' nanmax(matrixData(:,ii))];
            else
                [N,edges] = histcounts(matrixData(:,ii));
            end
            for jj=1:numel(N)
                
                % Get the subpopulation
                curIdx = matrixData(:,ii)>=edges(jj) & matrixData(:,ii)<edges(jj+1);
                % Calculate the width according to N(jj)
                xData(curIdx) = ii - N(jj)/Nmax*scatterWidth/2 + N(jj)/Nmax*scatterWidth*(rand(size(xData(curIdx),1),1));
            end
            if horizontalPlot
                scatter(ax, matrixData(:,ii), xData,'filled','MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','SizeData',markerSize)
            else
                scatter(ax, xData,matrixData(:,ii),'filled','MarkerFaceColor',[.3 .3 .3],'MarkerEdgeColor','none','SizeData',markerSize)
            end
        end
        hold on
    end
    alpha(.5)
end
if all(onlyOneDataPerEachGroup)
    onlyOneDataAllGroups=true;
end
if ~onlyOneDataAllGroups
    if notchOn %min(sum(~isnan(matrixData),1))>20 || 
        if horizontalPlot
            boxplot(ax,matrixData,'whisker',whiskerRatio,'notch','on',...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',1,'colors','k','Orientation','horizontal');%, 'labelorientation','inline');
        else
            boxplot(ax,matrixData,'whisker',whiskerRatio,'notch','on',...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',1,'colors','k');%, 'labelorientation','inline');
        end
    else % if the data is too small, don't use notch
%         boxplot(matrixData,'whisker',whiskerRatio*0.95,'notch','off',...
%             'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',0,'colors','k');%, 'labelorientation','inline');
        if horizontalPlot
            boxplot(ax,matrixData,...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',1,'colors','k','Orientation','horizontal');%, 'labelorientation','inline');
        else
            boxplot(ax,matrixData,...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',0,'colors','k');%, 'labelorientation','inline');
        end
    end
else
    xlim([0 numCategories+1])
    set(ax,'XTick',1:numCategories)
    set(ax,'XTicklabel',nameListNew)
end    
set(findobj(ax,'LineStyle','--'),'LineStyle','-')
set(findobj(ax,'tag','Median'),'LineWidth',2)
% set(gca,'XTick',1:numel(nameList))
% set(gca,'XTickLabel',nameList)
if ~horizontalPlot
    set(ax,'XTickLabelRotation',45)
end

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
                if (p<0.05 && forceShowP~=2) || forceShowP==1 
                    q=q+lineGap;
                    if horizontalPlot
                        line(ax,ones(1,2)*(maxPoint2+q),[k ii], 'Color','k')    
                    else
                        line(ax,[k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
                    end
                    q=q+lineGap;
                    if horizontalPlot
                        t=text(ax,maxPoint2+q, floor((k+ii)/2)+0.3, ['p=' num2str(p,'%3.2e') ' (r)']);
                        t.Rotation=45;
                    else
                        text(ax,floor((k+ii)/2)+0.3, maxPoint2+q,['p=' num2str(p,'%3.2e') ' (r)'])
                    end
                end
            else
                [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
                if (p<0.05 && forceShowP~=2) || forceShowP==1
                    q=q+lineGap;
                    if horizontalPlot
                        line(ax,ones(1,2)*(maxPoint2+q),[k ii], 'Color','k')    
                    else
                        line(ax,[k ii], ones(1,2)*(maxPoint2+q),'Color','k')    
                    end
                    q=q+lineGap;
                    if horizontalPlot
                        t=text(ax, maxPoint2+q, floor((k+ii)/2)+0.3, ['p=' num2str(p,'%3.2e') ' (t)']);
                        t.Rotation=45;
                    else
                        text(ax,floor((k+ii)/2)+0.3, maxPoint2+q,['p=' num2str(p,'%3.2e') ' (t)'])
                    end
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
    if horizontalPlot
        xlim([minPoint2-lineGap*10 maxPoint2+lineGap*10])
    else
        ylim([minPoint2-lineGap*10 maxPoint2+lineGap*10])
    end
else
    if horizontalPlot
        xlim auto
    else
        ylim auto
    end
end

set(ax,'FontSize',7)
set(findobj(ax,'Type','Text'),'FontSize',6)
sucess=true;
hold off

