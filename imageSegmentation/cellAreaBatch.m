% cellAreaBatch.m is a script that collects cell area from Segmentation.
% It works similarly with adhesionAnalysisBatch and strainEnergyBatch.
% Sangyoon Han Apr 2022
%% open necessary MLs
clear

[pathAnalysisAll, MLNames, groupNames, usedSelectedFoldersMat,...
    specificName,~,MLdirect]=chooseSelectedFolders;
% Asking user
disp('The current names are: ')
nameList = groupNames';
disp(nameList)
namesOK = input('Do you want to rename your condition names? (y/n)','s');
if strcmp(namesOK, 'y')
    for ii=1:numel(nameList)
        curName = input(['For ' nameList{ii} ': '], 's');
        if ~isempty(curName)
            nameList{ii} = curName;
        end
    end
    specificName = strjoin(nameList, '_');
end
%% Just in case when there is a larger condition.
largerCondition = input('Do you want to add a larger condition name, e.g., control or corona-treated etc (y/n)?','s');
if strcmp(largerCondition, 'y')
    largerConditionName = input('Enter the condition name:', 's');
    specificName = [largerConditionName specificName];
end

%% Output
rootAnalysis = fileparts(pathAnalysisAll{1});
% rootAnalysis = pathAnalysisAll{1};
summaryPath = [rootAnalysis '/CellAreaSummary' specificName];
ii=0;
while exist(summaryPath, 'dir')
    ii=ii+1;
    summaryPath = [rootAnalysis '/CellAreaSummary' specificName num2str(ii)];
end
figPath = [summaryPath '/Figs'];
mkdir(figPath)
dataPath = [summaryPath '/Data'];
mkdir(dataPath)
save([rootAnalysis filesep 'selectedFoldersForArea_' specificName '.mat'], 'rootAnalysis','pathAnalysisAll','MLNames','groupNames')
%% Load movieLists for each condition
numConditions = numel(pathAnalysisAll);
for k=1:numConditions
    MLAll(k) = MovieList.load([pathAnalysisAll{k} filesep MLNames{k}]);
end
%% Get area from segmentation pacakge
N=zeros(numConditions,1);
cellAreaGroup = cell(numConditions,1);

sampleMovie = MLAll(1).movies_{1};

for ii=1:numConditions
    N(ii) = numel(MLAll(ii).movies_);

    cellAreaCond = zeros(N(ii),1);
    curML = MLAll(ii);
    
    % Combine data per each condition (1,2,3,4 for 3.4, 18, 100, 100Y, respectively)
    p=0;
    for k=1:N(ii)
        curArea = getCellArea(curML.movies_{k});
        for m=1:numel(curArea)
            p=p+1;
            cellAreaCond(p) = curArea;
        end
    end
    cellAreaGroup{ii} = cellAreaCond;
end

%% Plotting cell area
try
    cellAreaCell = cellfun(@(x) cell2mat(cellfun(@(y) y(:),x,'unif',false)), cellAreaGroup,'unif',false);
catch
    cellAreaCell = cellfun(@(x) cell2mat(x),cellAreaGroup,'unif',false);
end
% cellAreaCell = cellfun(@(x) x(:), cellAreaCell,'unif',false);
h1=figure; 
boxPlotCellArray(cellAreaCell,nameList,1,1,1)
ylabel('Cell area (\mum^2)')
title('Cell area')
hgexport(h1,strcat(figPath,'/cellArea'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/cellArea'),'-v7.3')
tableCellarea=table(cellAreaCell,'RowNames',nameList);
writetable(tableCellarea,strcat(dataPath,'/CellArea.csv'))
%% error bar plot
errorBarPlotCellArray(cellAreaCell,nameList,1);
ylabel('Cell area (\mum^2)')
title('Cell area')
hgexport(h1,strcat(figPath,'/cellAreaScatter'),hgexport('factorystyle'),'Format','eps')
hgsave(h1,strcat(figPath,'/cellAreaScatter'),'-v7.3')

%% bar plot
