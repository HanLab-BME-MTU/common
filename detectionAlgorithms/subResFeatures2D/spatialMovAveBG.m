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

% mkitti: not necessary since we build this up
%allocate memory for output
% bgMean = NaN(imageSizeX,imageSizeY);
% bgStd = bgMean;

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
% bgMean(startPixelX-6+(1:size(bgMean2,1)),startPixelY-6+(1:size(bgMean2,2)) ) = bgMean2;
bgStdFull = reshape(bgStdFull,offsetSize);
bgStdFull = imresize(bgStdFull,blockSize,'nearest');
% bgStd(startPixelX-6+(1:size(bgStd2,1)),startPixelY-6+(1:size(bgStd2,2))) = bgStd2;
% 

% Conserve memory, what's the performance hit?
clear nhoodIdx

%% For loops eliminated by bsxfun magic, 4X optimization

%go over all pixels within limits
% for iPixelX = startPixelX : 11 : endPixelX
%     for iPixelY = startPixelY : 11 : endPixelY
%         
%         %get local image
%         imageLocal = imageLast5(iPixelX-15:min(iPixelX+15,imageSizeX),iPixelY-15:min(iPixelY+15,imageSizeY),:);
% %         imageLocal2 = imageLast5(nhoodIdx+(iPixelX-16)+(iPixelY-16)*imageSizeX);
%         %estimate robust mean and std
%         %first remove NaNs representing cropped regions
%         imageLocal = imageLocal(~isnan(imageLocal));
%         if ~isempty(imageLocal)
%             [bgMean1,bgStd1] = robustMean(imageLocal(:));
%             bgStd1 = max(bgStd1,eps);
%         else
%             bgMean1 = NaN;
%             bgStd1 = NaN;
%         end
%         
%         %put values in matrix representing image
%         bgMean(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5) = bgMean1;
%         bgStd(iPixelX-5:iPixelX+5,iPixelY-5:iPixelY+5) = bgStd1;
%         
%     end
% end

%% Finds are not necessary
%find limits of actual pixels filled up above
% firstFullX = find(~isnan(bgMean(:,startPixelY)),1,'first');
% lastFullX = find(~isnan(bgMean(:,startPixelY)),1,'last');
% firstFullY = find(~isnan(bgMean(startPixelX,:)),1,'first');
% lastFullY = find(~isnan(bgMean(startPixelX,:)),1,'last');

%% Map full values to original matrix size
firstFullX = startPixelX - blockRadius;
% lastFullX = iPixelX + 5;
lastFullX = firstFullX + size(bgMeanFull,1) - 1;
firstFullY = startPixelY - blockRadius;
% lastFullY = iPixelY + 5;
lastFullY = firstFullY + size(bgMeanFull,1) - 1;

%% Pad array by replication
% 
% %patch the rest, mkitti: Patch = Pad
% for iPixelY = firstFullY : lastFullY
%     bgMean(1:firstFullX-1,iPixelY) = bgMean(firstFullX,iPixelY);
%     bgMean(lastFullX+1:end,iPixelY) = bgMean(lastFullX,iPixelY);
%     bgStd(1:firstFullX-1,iPixelY) = bgStd(firstFullX,iPixelY);
%     bgStd(lastFullX+1:end,iPixelY) = bgStd(lastFullX,iPixelY);
% end
% for iPixelX = 1 : imageSizeX
%     bgMean(iPixelX,1:firstFullY-1) = bgMean(iPixelX,firstFullY);
%     bgMean(iPixelX,lastFullY+1:end) = bgMean(iPixelX,lastFullY);
%     bgStd(iPixelX,1:firstFullY-1) = bgStd(iPixelX,firstFullY);
%     bgStd(iPixelX,lastFullY+1:end) = bgStd(iPixelX,lastFullY);
% end
% Padding is asymmetric because only full blocks are considered
bgMean = padarray(bgMeanFull,[firstFullX firstFullY]-1,'replicate','pre');
bgMean = padarray(bgMean,[imageSizeX imageSizeY] - [lastFullX lastFullY],'replicate','post');
clear bgMeanFull

bgStd = padarray(bgStdFull,[firstFullX firstFullY]-1,'replicate','pre');
bgStd = padarray(bgStd,[imageSizeX imageSizeY] - [lastFullX lastFullY],'replicate','post');

end
