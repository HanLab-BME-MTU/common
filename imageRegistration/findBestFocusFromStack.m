function [averagingRange] = findBestFocusFromStack(curRefBeadStack,thresVariance,applySobel,method)
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
if nargin<4
    method = 'var';
end
% Per image
% midPixelsAllFrames = cell(refMD.zSize_,1);
if strcmp(method,'var')
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
elseif strcmp(method,'amp')
    % I will take 100th to 200th brightest pixels to guess the best
    % focus (usually the very brightest point is from one
    % extraordinary bead) - SH 20171010
    % Get 100th to 200th pixels
    midPixelsAllFrames = cell(numZ,1);
    for jj=1:numZ
        curRefImageFrame = curRefBeadStack(:,:,jj);
        curRefImageFrameSorted = sort(curRefImageFrame(:),'descend');
        midPixelsAllFrames{jj} = curRefImageFrameSorted(100:300);
    end
    meanMidInten = cellfun(@mean,midPixelsAllFrames);
    % Take top five frames
    [~,meanMidIntenIDs]=sort(meanMidInten,'descend');
    averagingRange = meanMidIntenIDs(1:5);
    % maxProf = reshape(max(max(curRefBeadStack)),[],1);
    % [~,maxIntenFrame]=max(maxProf);
    % minFocusedFrame=max(1,maxIntenFrame-2);
    % maxFocusedFrame=max(refMD.zSize_,maxIntenFrame+2);    
end

% This is now obsolete

