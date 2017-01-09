function [bgMean,bgStd] = spatialMovAveBG(imageLast5,imageSizeX,imageSizeY)

%the function in its current form assigns blocks of 11x11 pixels the
%same background values, for the sake of speed

% Khuloud Jaqaman
% Mark Kittisopikul, refactored and optimized October 2015

%TODO: Consider padding before neighborhood filtering

%% Setup values 
% constant parameters
blockSize = 11;
nhoodSize = 31;

blockRadius = (blockSize-1)/2; %5

%define pixel limits where moving average can be calculated
startPixelX = (nhoodSize+1)/2; %16
endPixelX = max(imageSizeX - startPixelX+1,startPixelX);
startPixelY = startPixelX;
endPixelY = max(imageSizeY - startPixelY+1,startPixelY);

%% Setup neighborhood indices
% Compute 2D neighborhood indices
nhoodIdx = bsxfun(@plus,(1:nhoodSize)',(0:nhoodSize-1)*imageSizeX);
% Extend 2D neighborhood to cover 3D
nhoodIdx = bsxfun(@plus,nhoodIdx,shiftdim((0:size(imageLast5,3)-1)*imageSizeX*imageSizeY,-1));
% Offsets for row and col translation of the neighborhood
offsets = bsxfun(@plus,(0:blockSize:endPixelX-startPixelY)',(0:blockSize:endPixelY-startPixelY)*imageSizeX);
% Sliding neighborhoods by column
nhoodIdx = bsxfun(@plus,nhoodIdx(:),offsets(:)');
offsetSize = size(offsets);
clear offsets

%% Get robust mean and std for each neighborhood
% Calculate robust mean and std
[bgMeanFull,bgStdFull] = robustMean(imageLast5(nhoodIdx));
bgMeanFull = reshape(bgMeanFull,offsetSize);
bgMeanFull = imresize(bgMeanFull,blockSize,'nearest');
bgStdFull = reshape(bgStdFull,offsetSize);
bgStdFull = imresize(bgStdFull,blockSize,'nearest');

% Conserve memory, what's the performance hit?
clear nhoodIdx

%% Map full values to original matrix size
firstFullX = startPixelX - blockRadius;
lastFullX = firstFullX + size(bgMeanFull,1) - 1;
firstFullY = startPixelY - blockRadius;
lastFullY = firstFullY + size(bgMeanFull,2) - 1;

%% Pad array by replication
% 
% %patch the rest, mkitti: Patch = Pad
% Padding is asymmetric because only full blocks are considered
bgMean = padarray(bgMeanFull,[firstFullX firstFullY]-1,'replicate','pre');
bgMean = padarray(bgMean,[imageSizeX imageSizeY] - [lastFullX lastFullY],'replicate','post');
clear bgMeanFull

bgStd = padarray(bgStdFull,[firstFullX firstFullY]-1,'replicate','pre');
bgStd = padarray(bgStd,[imageSizeX imageSizeY] - [lastFullX lastFullY],'replicate','post');

end
