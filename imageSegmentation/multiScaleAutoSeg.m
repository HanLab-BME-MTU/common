function multiScaleAutoSeg(MD, iChan, varargin)
% MSA_Seg (multi-scale automatic segmentation) 
% Segment a single cell image by combining segmentation
% results obtained at multiple smoothing scales. Since it requires only one
% tuning parameters (tightness) and ‘tightness’=0.5 works well for many cases, 
% it achieves almost automatic segmentation.
% Usage :
%           MSA_Seg(MD, 1, 'tightness', 0.5, 'imagesOut', 1)
% 2017/05/29, Jungsik Noh
% Updated Andrew R. Jamieson Sept. 2017


ip = inputParser;
ip.addParameter('type', 'middle', @(x) strcmp('middle',x) || strcmp('tight', x) || strcmp('minmax', x));
ip.addParameter('tightness', -1);
ip.addParameter('imagesOut', 0);
% ip.addParameter('parpoolNum', 1);

ip.parse(varargin{:});
p = ip.Results;

tic

%% -------- Parameters ---------- %%
dName = ['MSASeg_refined_masks_channel_', num2str(iChan)];
% String for naming the mask directories for each channel

outDir = fullfile(MD.outputDirectory_, 'MultiScaleAutoSeg', dName);
if ~isdir(outDir); mkdir(outDir); end
    
pString = ['MSA_', 'refined_'];      %Prefix for saving masks to file

imFileNamesF = MD.getImageFileNames(iChan);
imFileNames = imFileNamesF{1};


%% Load images
I = cell(MD.nFrames_, 1);
imgStack = [];
for fr=1:MD.nFrames_
    I{fr} = MD.channels_(iChan).loadImage(fr);
    imgStack = cat(3, imgStack, I{fr});
end

%%% Segmentation Diagnostics
imageStats(imgStack, MD);

%% Multi Scale Segmentation
refinedMask = cell(MD.nFrames_, 1);
currType = p.type;
currTightness = p.tightness;

parfor fr = 1:MD.nFrames_
    disp('=====');
    disp(['Frame: ', num2str(fr)]);    
    im = I{fr};
    refinedMask{fr} = multiscaleSeg_im(im, 'type', currType, 'tightness', currTightness);

    %Write the refined mask to file
    imwrite(refinedMask{fr}, fullfile(outDir, [pString, imFileNames{fr}]));
end

%% imagesOut
if p.imagesOut == 1
    
    if p.tightness >= 0
        prefname = ['tightness_', num2str(p.tightness)];
    else 
        prefname = p.type;
    end
       
    dName2 = ['MSASeg_maskedImages_', prefname,  '_channel_', num2str(iChan)];
    imOutDir = fullfile(MD.outputDirectory_, 'MultiScaleAutoSeg', dName2);
    if ~isdir(imOutDir); mkdir(imOutDir); end

    imgStack = [];
    for fr = 1:MD.nFrames_
        imgStack = cat(3, imgStack, I{fr});
    end    
    allint = imgStack(:);
    intmin = quantile(allint, 0.002);
    intmax = quantile(allint, 0.998);

    ftmp = figure;
    for fr = 1:MD.nFrames_
        figure(ftmp);
        imshow(I{fr}, [intmin, intmax]);
        hold on;
        bdd = bwboundaries(refinedMask{fr});
        bdd1 = bdd{1};
        plot(bdd1(:,2), bdd1(:,1), 'r');
        hold off;
        
        h = getframe(gcf);
        imwrite(h.cdata, fullfile(imOutDir, imFileNames{fr}), 'tif')
    end
end

%%
toc
disp('Multi-Scale Automatic Segmentation is done!')

end


function imageStats(imgStack, MD)
    %%% Segmentation Diagnostics
    %% TS of 5 numbers
    pixelmat = reshape(imgStack, [], MD.nFrames_);
    pixelmat1 = pixelmat;
    pixelmat1(pixelmat1 == 0) = NaN;

    mts = mean(pixelmat1, 1, 'omitnan');
    medts = median(pixelmat1, 1, 'omitnan');
    q1ts = quantile(pixelmat1, 0.25, 1);
    q3ts = quantile(pixelmat1, 0.75, 1);
    q99ts = quantile(pixelmat1, 0.99, 1);
    q01ts = quantile(pixelmat1, 0.01, 1);

    fts = figure;
    plot(mts)
    hold on

    plot(medts)
    plot(q1ts)
    plot(q3ts)
    plot(q01ts)
    plot(q99ts)
    hold off

    legend('Mean', 'Median', 'Perct25', 'Perct75', 'Perct01', 'Perct99');
    title('Time series of 5 summary statistics');

    %%
    saveas(fts, fullfile(MD.outputDirectory_, 'MultiScaleAutoSeg', 'TS_of_5statistics.png'), 'png')
    saveas(fts, fullfile(MD.outputDirectory_, 'MultiScaleAutoSeg', 'TS_of_5statistics.fig'), 'fig')
end