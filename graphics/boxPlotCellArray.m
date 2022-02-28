function [sucess]=boxPlotCellArray(cellArrayData,nameList,convertFactor,notchOn,plotIndivPoint,varargin)
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
% Example:
% sampleArray{1}=randn(100,1);
% sampleArray{2}=randn(100,1)*1.3+2;
% sampleArray{3}=randn(50,1)*1.5+1;
% nameList={'G1','G2','G3'};
% boxPlotCellArray(sampleArray,nameList,1,1,1)
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
addParameter(ip,'forceShowP',false);
addParameter(ip,'markerSize',2);
addParameter(ip,'ax',gca);
addParameter(ip,'horizontalPlot',false);
addParameter(ip,'forceTtest',false);
addParameter(ip,'singleColor',0);
parse(ip,cellArrayData,nameList,convertFactor,notchOn,plotIndivPoint,varargin{:});
ax=ip.Results.ax;
horizontalPlot=ip.Results.horizontalPlot;
markerSize = ip.Results.markerSize;
forceShowP = ip.Results.forceShowP;
notchOn = ip.Results.notchOn;
plotIndivPoint = ip.Results.plotIndivPoint;
convertFactor = ip.Results.convertFactor;
forceTtest =  ip.Results.forceTtest;
singleColor = ip.Results.singleColor;

nameList(idEmptyData)=[];

boxWidth=0.5;
scatterWidth=boxWidth*2;
whiskerRatio=0.5;
matrixData=matrixData*convertFactor;
try
    nameListNew = cellfun(@(x,y) [x '(N=' num2str(sum(~isnan(y))) ')'],nameList,cellArrayData,'UniformOutput', false);
catch
    nameList=arrayfun(@num2str,(1:numel(cellArrayData))','unif',false);
    nameListNew = cellfun(@(x,y) [x '(N=' num2str(sum(~isnan(y))) ')'],nameList,cellArrayData,'UniformOutput', false);
end

% Generating color group
%      r  red       1               
%      g  green     2
%      b  blue      3
%      y  yellow    4            
%      m  magenta   5           
%      c  cyan      6        
%      a  apple green 7  
%      d  dark gray 8   
%      e  evergreen 9    
%      f  fuchsia   10                  
%      h  honey     11     
%      i  indigo    12                  
%      j  jade      13               
%      l  lilac     14            
%      n  nutbrown  15         
%      p  pink      16                  
%      q  kumquat   17             
%      s  sky blue  18                    
%      t  tan       19                          
%      u  umber     20                            
%      v  violet    21                           
%      z  zinc      22
%      k  black     23
colorSwitch = 'qkjlrabnfilkpuvdt';
color = extendedColors(colorSwitch);

onlyOneDataAllGroups=false;
numCategories = numel(cellArrayData);
while size(color,1)<numCategories
    color = [color; color];
end

onlyOneDataPerEachGroup=false(1,numCategories);
if plotIndivPoint
    % individual data jitter plot
%     numCategories = size(matrixData,2);
    % width=boxWidth/2;
    uArrayAll = cell(size(cellArrayData));
    Nall = cell(size(cellArrayData));
    for ii=1:numCategories
        Nall{ii} = histcounts(matrixData(:,ii));
        uArrayAll{ii} = unique(matrixData(~isnan(matrixData(:,ii)),ii));
    end
%     if mean(cellfun(@(x,y) length(y)/length(x),uArrayAll,cellArrayData))>3
%         Nmax = max(cellfun(@max,cellfun(@(x,y) arrayfun(@(z) sum(y==z),x) ,uArrayAll,cellArrayData,'unif',false)));
%     else
        Nmax = max(cellfun(@max,Nall));
%     end
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
                curMaxNratio = max(N/Nmax);
                scatterWidth = boxWidth/curMaxNratio;
            else
                [N,edges] = histcounts(matrixData(:,ii));
            end
            for jj=1:numel(N)
                
                % Get the subpopulation
                curIdx = matrixData(:,ii)>=edges(jj) & matrixData(:,ii)<edges(jj+1);
                % Calculate the width according to N(jj)
                xData(curIdx) = ii - N(jj)/Nmax*scatterWidth/2 + N(jj)/Nmax*scatterWidth*(rand(size(xData(curIdx),1),1));
            end
            if singleColor==0
                curColor = color(ii,:);
            else
                curColor = singleColor;
            end
            if horizontalPlot
                scatter(ax, matrixData(:,ii), xData,'filled','MarkerFaceColor',curColor*0.5,'MarkerEdgeColor','none','SizeData',markerSize)
            else
                scatter(ax, xData,matrixData(:,ii),'filled','MarkerFaceColor',curColor*0.5,'MarkerEdgeColor','none','SizeData',markerSize)
            end
        end
        hold on
    end
    alpha(.5)
end
if all(onlyOneDataPerEachGroup)
    onlyOneDataAllGroups=true;
end

% color setting for box plots
if singleColor==0
    colors = color(ii,:);
else
    colors = singleColor;
end

if ~onlyOneDataAllGroups %Here I need to change color option for single color.
    if notchOn %min(sum(~isnan(matrixData),1))>20 || 
        if horizontalPlot
            boxplot(ax,matrixData,'whisker',whiskerRatio,'notch','on',...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',0,'colors',colors,'Orientation','horizontal');%, 'labelorientation','inline');
        else
            boxplot(ax,matrixData,'whisker',whiskerRatio,'notch','on',...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',0,'colors',colors);%, 'labelorientation','inline');
        end
    else % if the data is too small, don't use notch
%         boxplot(matrixData,'whisker',whiskerRatio*0.95,'notch','off',...
%             'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',0,'colors','k');%, 'labelorientation','inline');
        if horizontalPlot
            boxplot(ax,matrixData,'whisker',whiskerRatio,...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',1,'colors',colors,'Orientation','horizontal');%, 'labelorientation','inline');
        else
            boxplot(ax,matrixData,'whisker',whiskerRatio,...
                'labels',nameListNew,'symbol','','widths',boxWidth,'jitter',0,'colors',colors);%, 'labelorientation','inline');
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
% We assume that at each cell array, there is Nx1 vectors (not 1xN vectors)
% we flip them if they are in row vectors:
if size(cellArrayData{1},1)==1 && size(cellArrayData{1},2)>1
    cellArrayData = cellfun(@transpose, cellArrayData, 'unif', false);
elseif numel(cellArrayData)>1 && size(cellArrayData{2},1)==1 && size(cellArrayData{2},2)>1
    cellArrayData = cellfun(@transpose, cellArrayData, 'unif', false);
end
% We also have to format cellArrayData to be only a colume.
if size(cellArrayData,1)==1 && size(cellArrayData,2)>1
    cellArrayData = cellArrayData';
end
% hold on
% perform ranksum test for every single combination
maxPoint = quantile(matrixData,[0.25 0.75]);
if size(matrixData,2)>1 && size(matrixData,1)>1
    maxPoint2 = maxPoint(2,:)+(maxPoint(2,:)-maxPoint(1,:))*whiskerRatio;
    maxPoint2 = max(maxPoint2);
else
    maxPoint2 = maxPoint(2)+(maxPoint(2)-maxPoint(1))*whiskerRatio;
end
lineGap=maxPoint2*0.05;
xGap = 0.05;
% q=-2*lineGap;
method='';
qUsed=[];
for k=1:(numConditions-1)
%     q=-2*lineGap + 0.5*lineGap*(k-1); %0.5*lineGap*mod(k-1,2);
    for ii=k+1:numConditions
        if numel(cellArrayData{k})>1 && numel(cellArrayData{ii})>1
            if kstest(cellArrayData{k}) && ~forceTtest % this means the test rejects the null hypothesis
                method = 'ranksum test';
                [p]=ranksum(cellArrayData{k},cellArrayData{ii});
            else
                method = 'unpaired t-test';
                [~,p]=ttest2(cellArrayData{k},cellArrayData{ii});
            end
            if (p<0.05 && forceShowP~=2) || forceShowP==1 
%                 q = quantile(cell2mat(cellArrayData(k:ii)),0.92);
                q = max(cellfun(@(x) quantile(x,0.95),cellArrayData(k:ii)));
                tol=0.95*lineGap;

                a=find(abs(qUsed-q(1))<tol,1);
                while ~isempty(a) %ismember(q,qUsed) 
                    q=qUsed(a)+lineGap;
                    a=find(abs(qUsed-q)<tol,1);
                end
                qUsed = [qUsed q];
                q=q+lineGap;
                if horizontalPlot
                    line(ax,ones(1,2)*(maxPoint2+q),[k ii], 'Color','k')    
                else
                    line(ax,[k+xGap ii-xGap], ones(1,2)*(q),'Color','k')    
                end
                q=q+lineGap*0.5;
                if horizontalPlot
                    t=text(ax,maxPoint2+q, floor((k+ii)/2)+0.3, ['p=' num2str(p,'%2.2e') ' (r)']);
                    t.Rotation=45;
                else
%                         text(ax,floor((k+ii)/2)+0.3, maxPoint2+q,['p=' num2str(p,'%2.2e') ' (r)'])
                    text(ax,((k+ii)/2)-0.2, q,['p=' num2str(p,2)]); %'%2.2e'
                end
            end
        end
    end
end
% q=q+lineGap*3;
% minPoint = quantile(matrixData,[0.25 0.75]);
% if size(matrixData,2)>1 && size(minPoint,1)>1
%     minPoint2 = minPoint(1,:)-(minPoint(2,:)-minPoint(1,:))*whiskerRatio;
%     minPoint2 = min(minPoint2);
% else
%     minPoint2 = minPoint(2)-(minPoint(2)-minPoint(1))*whiskerRatio;
% end
% if plotIndivPoint && ~forceShowP
%     maxPoint2 = quantile(matrixData(:),0.99);
% else
%     maxPoint2 = maxPoint2+q;
% end
% if maxPoint2>minPoint2-lineGap*2
%     if horizontalPlot
%         xlim([minPoint2-lineGap*10 maxPoint2+lineGap*10])
%     else
%         try
%             ylim([minPoint2-lineGap*5 maxPoint2+lineGap*12])
%         catch
%             ylim auto
%         end
%     end
% else
%     if horizontalPlot
%         xlim auto
%     else
%         ylim auto
%     end
% end

xlim auto
ylim auto

set(ax,'FontSize',7)
set(findobj(ax,'Type','Text'),'FontSize',6)
sucess=true;
text(ax,ii-1.5,ax.YLim(1)+lineGap,method)

hold off

