function compMatrices = discriminationMatrix(dataStructure)
%DISCRIMINATIONMATRIX generates discriminationMatrices comparing means and distributions of lists of values in the input data structures
%
% SYNOPSIS compMatrices = discriminationMatrix(dataStructure)
%
% INPUT    dataStructure : n-by-1 structure containing the lists of values
%                          that are to be compared
%                          Fieldnames can be arbitrary
%
% OUTPUT   compMatrices  : Structure with fieldnames equal to the
%                     fieldnames of the data structure, containing
%                     n-by-n matrices with p-values for
%                     - below the diagonal: T-test for equality of means
%                     - above the diagonal: Kolmogorov-Smirnov test for
%                         equality of the distributions, shifted by their
%                         means.
%                     Values in the diagonal are 1.01
%
% c: jonas, 12/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=====================
% TEST INPUT
%=====================

% nargin
if nargin == 0 || isempty(dataStructure) || ~isstruct(dataStructure)
    error('Please input a non-empty data structure in DISCRIMINATIONMATRIX')
end

% length structure
nGroups = length(dataStructure);
if nGroups < 2
    error('Please try to compare at least two sets of data!')
end

% number of data sets to compare
dataNames = fieldnames(dataStructure);
nData = length(dataNames);

% loop through all the data and get means (check whether empty, too)
groupMeans = zeros(nGroups,nData);
for iData = 1:nData
    for iGroup = 1:nGroups
        if ~isempty(dataStructure(iGroup).(dataNames{iData}))
            
            groupMeans(iGroup,iData) = ...
                mean(dataStructure(iGroup).(dataNames{iData}));
        else
            error('no empty entries are allowed in the dataStructure')
        end
    end
end

%======================



%=====================
% TEST DATA
%=====================

% preassign output
rawMat = repmat(1.01,[nGroups,nGroups]);

% fore each data entry: Compare all the groups
for iData = 1:nData

    % preassign output
    compMatrices.(dataNames{iData}) = rawMat;

    for gi = 1:nGroups
        for gj = gi-1:-1:1
            % compare mean:  ttest. Below diagonal
            [dummy,compMatrices.(dataNames{iData})(gi,gj)] = ...
                ttest2(...
                dataStructure(gi).(dataNames{iData}),...
                dataStructure(gj).(dataNames{iData}),...
                0.05,'both','unequal');
            % compare distribution (shifted data): kolmogorov-smirnov test.
            % Above diagonal
            [dummy,compMatrices.(dataNames{iData})(gj,gi)] = ...
                kstest2(...
                dataStructure(gi).(dataNames{iData})-groupMeans(gi,iData),...
                dataStructure(gj).(dataNames{iData})-groupMeans(gj,iData));
        end
    end
end
