function [averagingRange] = findBestFocusFromStack(curRefBeadStack,thresVariance,applySobel)
% function [averagingRange] =
% findBestFocusFromStack(curRefBeadStack,thresVariance,applySobel)
% identifies the best focus by top most variance among the stack images.
% Output is indices for the stack images. Once applySobel is true, we apply
% the sobel filter to the images and get the variance. The output is
% determined by whatever whose varance, normalized by max variance, is over
% thresariance.
% Sangyoon Han, Sep 2018

numZ = size(curRefBeadStack,3);
curVar = zeros(numZ,1);
% Per image
% midPixelsAllFrames = cell(refMD.zSize_,1);
for jj=1:numZ
    curRefImageFrame = curRefBeadStack(:,:,jj);
    % Apply sobel filter
    if applySobel
        H = fspecial('sobel');
        imgSobeled = imfilter(curRefImageFrame,H,'replicate'); 
    end
    % Get the variance
    if applySobel
        curVar(jj) = var(double(imgSobeled(:)));
    else
        curVar(jj) = var(double(curRefImageFrame(:)));
    end
end
maxVar = max(curVar);
normVar = curVar/maxVar;
averagingRange = find(normVar>=thresVariance);

end

% This is now obsolete
% % I will take 100th to 200th brightest pixels to guess the best
% % focus (usually the very brightest point is from one
% % extraordinary bead) - SH 20171010
% % Get 100th to 200th pixels
% midPixelsAllFrames = cell(refMD.zSize_,1);
% for jj=1:refMD.zSize_
%     curRefImageFrame = curRefBeadStack(:,:,jj);
%     curRefImageFrameSorted = sort(curRefImageFrame(:),'descend');
%     midPixelsAllFrames{jj} = curRefImageFrameSorted(100:300);
% end
% meanMidInten = cellfun(@mean,midPixelsAllFrames);
% % Take top five frames
% [~,meanMidIntenIDs]=sort(meanMidInten,'descend');
% averagingRange = meanMidIntenIDs(1:5);
% % maxProf = reshape(max(max(curRefBeadStack)),[],1);
% % [~,maxIntenFrame]=max(maxProf);
% % minFocusedFrame=max(1,maxIntenFrame-2);
% % maxFocusedFrame=max(refMD.zSize_,maxIntenFrame+2);
