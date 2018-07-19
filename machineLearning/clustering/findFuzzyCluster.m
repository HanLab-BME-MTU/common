function index = findFuzzyCluster(data,nCluster)
%This is a wrapper for the fcm function. It parses the output.
%
%Usage: 
%       index = findFuzzyCluster(data,nCluster)
%
%Input:
%       data     - matrix(#ofFeatures,#ofPoints)
%       nCluster - number of clusters
%
%Output:
%       index - indexes of the input for each cluster
%
%Marco Vilela, 2013

[~,U] = fcm(data, nCluster);
maxU  = max(U);
index = cell(1,nCluster);
meanV = nan(1,nCluster);

for iC = 1:nCluster
    index{iC} = find(U(iC,:) == maxU);
    meanV(iC) = mean(data(index{iC}));
end

[~,idx] = sort(meanV);
index   = index(idx);

end
