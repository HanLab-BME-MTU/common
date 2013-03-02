function [fracClustTotalMS,fracClustPerFrameMS,fracClustTotalMovie,fracClustPerFrameMovie] = calcFracClust(tracksAggreg,maxSize)
%CALCFRACCLUST calculates the fraction of clusters of various sizes from series of tracks
%
%SYNOPSIS [fracClustPerFrameMovie,fracClustTotalMovie] = calcFracClust(tracksAggreg,masSize)
%
%INPUT  tracksAggreg : Cell array of output of aggregStateFromCompTracks.
%                      1 entry = 1 movie/simulation.
%       maxSize      : Maximum cluster size to look into. Any clusters
%                      larger than maxSize will be lumped into the bin for
%                      maxSize.
%                      Optional. Default: maximum cluster size in data.
%
%OUTPUT 
%
%Khuloud Jaqaman, February 2013

%% Input

%get number of movies
numMovie = length(tracksAggreg);

%get matrices of aggregation state
aggregMat = cell(numMovie,1);
for iMovie = 1 : numMovie
    [~,~,~,~,aggregMat{iMovie}] = convStruct2MatIgnoreMS(tracksAggreg{iMovie}.defaultFormatTracks);
end

%get number of frames
aggregMatAll = vertcat(aggregMat{:});
numFrames = size(aggregMatAll,2);

%get maximum cluster size if not input
if nargin < 2 || isempty(maxSize)
    maxSize = max(aggregMatAll(:));
end

%% Cluster fractions

%fractions per frame in each movie
fracClustPerFrameMovie = cell(numMovie,1);
for iMovie = 1 : numMovie
    fracClust = hist(aggregMat{iMovie},1:maxSize);
    fracClust = fracClust ./ repmat(sum(fracClust),maxSize,1);
    fracClustPerFrameMovie{iMovie} = fracClust;
end

%overall fractions in each movie
fracClustTotalMovie = cell(numMovie,1);
for iMovie = 1 : numMovie
    fracClust = hist(aggregMat{iMovie}(:),1:maxSize);
    fracClust = fracClust / sum(fracClust);
    fracClustTotalMovie{iMovie} = fracClust;
end

%mean and std of fractions per frame
fracClustPerFrameAll = vertcat(fracClustPerFrameMovie{:});
fracClustPerFrameMS = zeros(maxSize,numFrames,2);
for iSize = 1 : maxSize
    fracClustPerFrameMS(iSize,:,1) = mean(fracClustPerFrameAll(iSize:maxSize:end,:));
    fracClustPerFrameMS(iSize,:,2) = std(fracClustPerFrameAll(iSize:maxSize:end,:));
end
    
%mean and std of overall fractions
fracClustTotalAll = vertcat(fracClustTotalMovie{:});
fracClustTotalMS = [mean(fracClustTotalAll); std(fracClustTotalAll)]';



%% ~~~ the end ~~~
