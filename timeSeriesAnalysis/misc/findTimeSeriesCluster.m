function out = findTimeSeriesCluster(TS,fVector,clusterMethod,varargin)
%This function cluster time series data using a set of different features.
%
%USAGE:
%       out = findTimeSeriesCluster(TS,fVector,clusterMethod,varargin)
%
%Input:
%       TS      - cell array where each element is a TS matrix
%       fVector - cell array where each element is a vector with the features ID
%                  Features ID:
%                                1 - Time series energy
%                                2 - Time series mean
%                                3 - Time series variance
%                                4 - Time series ACF
%                                5 - Time series CCF
%
%       clusterMethod - scalar ID for the cluster method
%                                1 - K-means 
%                                2 - fuzzy K-means
%                                3 - Mean shift
%
%       Optional:
%                maxLag     - ACF or CCF maximum lag, if 4 and/or 5 is chosen in fVector
%                clusterSet - cell array with the input options for each clustering method. Help the cluster method function to figure this out out 
%
%Output:
%       out - cell array where each element contains the indexes for each cluster
%
%Marco Vilela, 2013


%% Parsing input

ip = inputParser;
ip.addRequired('TS',@(x) iscell(x) || ismatrix(x));
ip.addRequired('fVector',@isvector );
ip.addRequired('clusterMethod',@isscalar);

ip.addParamValue('maxLag',   0,@isscalar);
ip.addParamValue('clusterSet', {2,'Distance','sqEuclidean','Replicates',10}, @iscell);                 

ip.parse(TS,fVector,clusterMethod,varargin{:});
maxLag     = ip.Results.maxLag;
clusterSet = ip.Results.clusterSet;
featInput  = {'maxLag',maxLag};

if numel(TS) ~= numel(fVector)
    error('Number of Variables and number of feature vector does not match')
end


%% Clustering methods

methods{1} = @(x,y) kmeans(x,y{:});
methods{2} = @(x,y) findFuzzyCluster(x,y{:});
methods{3} = @(x,y) MeanShiftClustering(x,y{:});

%% estimating the features

featMatrix = cell2mat( cellfun(@(x,y) gettingFeatures(x,y,featInput),TS,fVector,'Unif',0)' );

%% Running clustering algorithm

out = methods{clusterMethod}(featMatrix',clusterSet);

%% Formatting output

if clusterMethod == 1
    
    out = cellfun(@(x) find(out == x),num2cell(1:clusterSet{1}),'Unif',0);
    
elseif clusterMethod == 3
    
    out = arrayfun(@(x) x.ptIdData,out,'Unif',0);
end

end

function out = gettingFeatures(TS,fVector,maxLag)

% Time Series Features

feature{1} = @(x)   x*x';%Energy
feature{2} = @(x)   nanmean(x);%Mean
feature{3} = @(x)   nanvar(x);%Variance
feature{4} = @(x)   nanCrossCorrelation(x,x,'maxLag',maxLag);%ACF
feature{5} = @(x,y) nanCrossCorrelation(x,y,'maxLag',maxLag);%CC

selectedFeat = @(x) cellfun(@(y) y(x),feature(fVector),'Unif',0);
estFeat      = cellfun(@(x) selectedFeat(x),num2cell(TS,2),'Unif',0);
out          = cell2mat(cellfun(@(x) cell2mat(x'),estFeat,'Unif',0)');

end
