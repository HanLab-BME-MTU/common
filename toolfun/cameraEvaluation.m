function [averageImage,sdImage,hotImage,xPattern, yPattern, xPattern2, yPattern2]=cameraEvaluation(darkFrames)
%CAMERAEVALUATION tests cameras for static and dynamic patterns
%
% The code performs an analysis on a series of dark images. 
%       First, the average and standard deviation of the image series are
%       calculated. These are the static patterns. For better visibility,
%       hot pixels in the static patterns are masked with the robust mean
%       in the images.
%       To find dynamic patterns, the dark images are corrected by the
%       average and divided by the standard deviation. Then, every single
%       image is analyzed twice for dynamic correlation: Once by looking
%       for correlation on the entire corrected images, and once on images
%       with a 20% border removed. For every lag, the number of significant
%       correlations is counted (95% level). If there are more than 5%
%       significant lags in any given direction (=above the red line in the
%       graphs), there is a correlation.
%       The code writes the number of hot pixels to the commandline.
%
% SYNOPSIS testCamera(darkFrames)
%
% INPUT    darkFrames : 3D, 4D or 5D image stack (will be converted to a 3D
%                       image stack) 
%
% OUTPUT   averageImage : average of the image series
%          sdImage      : standard deviation of the image series
%          hotImage     : binary image of hot pixels
%          xPattern, yPattern, xPattern2, yPattern2
%                       : Logical arrays with the incidence of significant
%                         correlations.
%
% c: jonas, 10/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters
borderPercent = 20;

% convert image
movieSize = size(darkFrames);
darkFrames = reshape(darkFrames,movieSize(1),movieSize(2),[]);
movieSize = size(darkFrames);

% take average and std
averageImage = mean(darkFrames,3);
sdImage = std(darkFrames,0,3);

[robMean, dummy, inlierIdx] = robustMean(averageImage(:));
[robStd] = robustMean(sdImage(:));
outlierIdx = missingIndices(inlierIdx,prod(movieSize(1:2)));
hotImage = zeros(movieSize(1:2));
hotImage(outlierIdx) = 1;
hotImage = logical(hotImage);
numberOfHotPixels = nnz(hotImage);
ratioOfHotPixels = numberOfHotPixels/prod(movieSize(1:2));


% correct averageImage, stdImage for hot pixels
averageImageC = averageImage;
averageImageC(hotImage) = robMean;
sdImageC = sdImage;
sdImageC(hotImage) = robStd;

uiViewPanel;
imshow(averageImageC,[]);
set(gcf,'Name','Average background');
colormap('jet')


uiViewPanel;
imshow(sdImageC,[]);
set(gcf,'Name','Std background');
colormap('jet')

figure('Name','Average background (no hot pixels)'),
histogram(averageImageC(~hotImage));
figure('Name','Std background (no hot pixels)');
histogram(sdImageC(~hotImage));

% subtract static pattern, divide by sd
darkFrames = darkFrames - repmat(averageImage,[1,1,movieSize(3)]);
darkFrames = darkFrames ./ repmat(sdImage,[1,1,movieSize(3)]);

% plot average corrected darkFrame
% correctedAverageImage = mean(darkFrames,3);
% correctedSdImage = std(darkFrames,0,3);
% uiViewPanel,imshow(correctedAverageImage,[])
% set(gcf,'Name','average of corrected background')
% colormap('jet')
% uiViewPanel,imshow(correctedSdImage,[])
% set(gcf,'Name','std of corrected background')
% colormap('jet')
% 
% figure('Name','Average corrected background')
% histogram(correctedAverageImage(:));
% 
% figure('Name','Std corrected background')
% histogram(correctedSdImage(:));


% get dynamic pattern. Do either 10 times the number of pixels along the
% longer side of the image, or half the total number of pixels in the image
nCorr = min(max(10*movieSize(1:2)),ceil(prod(movieSize(1:2))/2));
[xPattern, yPattern] = deal(zeros(nCorr,movieSize(3)));
tic
for i = 1:movieSize(3)
    [xPattern(:,i), yPattern(:,i)] = ...
        normACfunc(darkFrames(:,:,i), -1);
    disp(sprintf('iteration: %i/%i, time: %9.3f',i,movieSize(3),toc))
end

% check for significance: is autocorrelation larger than 
% 1.96/sqrt(numberOfPixels)?
numberOfPixels = prod(movieSize(1:2));
xPattern = abs(xPattern) > 1.96/sqrt(numberOfPixels);
yPattern = abs(yPattern) > 1.96/sqrt(numberOfPixels);

% the first row is obviously not informative
xPattern(1,:) = logical(0);
yPattern(1,:) = logical(0);

% uiViewPanel,imshow(xPattern',[])
% uiViewPanel,imshow(yPattern',[])

% sum the incidence of significant lags. Normalize by the number of frames
sx = sum(xPattern,2)/movieSize(3);
sy = sum(yPattern,2)/movieSize(3);

figure('Name','Correlation X'),stem([0:nCorr-1]',sx);
hold on
plot([0:nCorr-1],0.05*ones(nCorr,1),'r')
figure('Name','Correlation Y'),stem([0:nCorr-1],sy);
hold on
plot([0:nCorr-1],0.05*ones(nCorr,1),'r')

% do again, but without 20% pixels all around
xPercent = floor(movieSize(1:2)/(borderPercent/100));
if all(xPercent > 2)
    
    darkFrames = darkFrames(xPercent(1)+1:end-xPercent(1),...
        xPercent(2)+1:end-xPercent(2),:);
    movieSize = size(darkFrames);
    nCorr = min(max(10*movieSize(1:2)),ceil(prod(movieSize(1:2))/2));
    [xPattern2, yPattern2] = deal(zeros(nCorr,movieSize(3)));
    tic
    for i = 1:movieSize(3)
        [xPattern2(:,i), yPattern2(:,i)] = ...
            normACfunc(darkFrames(:,:,i), -1);
        disp(sprintf('iteration: %i/%i, time: %9.3f',i,movieSize(3),toc))
    end

    % check for significance: is autocorrelation larger than
    % 1.96/sqrt(numberOfPixels)?
    numberOfPixels = prod(movieSize(1:2));
    xPattern2 = abs(xPattern2) > 1.96/sqrt(numberOfPixels);
    yPattern2 = abs(yPattern2) > 1.96/sqrt(numberOfPixels);

    % the first row is obviously not informative
    xPattern2(1,:) = logical(0);
    yPattern2(1,:) = logical(0);

    % uiViewPanel,imshow(xPattern',[])
    % uiViewPanel,imshow(yPattern',[])

    % sum the incidence of significant lags. Normalize by the number of frames
    sx = sum(xPattern2,2)/movieSize(3);
    sy = sum(yPattern2,2)/movieSize(3);

    figure('Name',...
        sprintf('Correlation X (minus 2x %i border pix)',xPercent(1)));
    stem([0:nCorr-1]',sx);
    hold on
    plot([0:nCorr-1],0.05*ones(nCorr,1),'r')
    figure('Name',...
        sprintf('Correlation Y (minus 2x %i border pix)',xPercent(2)));
    stem([0:nCorr-1],sy);
    hold on
    plot([0:nCorr-1],0.05*ones(nCorr,1),'r')
end


% display number of hot pixels
disp(sprintf('Number of hot pixels: %i (%2.3f)',numberOfHotPixels, ...
    100*ratioOfHotPixels));

