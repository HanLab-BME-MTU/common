function movieData = image_flatten(movieData, varargin)

% Find the process of segmentation mask refinement.
nProcesses = length(movieData.processes_);

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess==0
    msg('Please set parameters for Image Flatten.')
    return;
end

funParams=movieData.processes_{indexFlattenProcess}.funParams_;

selected_channels = funParams.ChannelIndex;
flatten_method_ind = funParams.method_ind;
Gaussian_sigma = funParams.GaussFilterSigma;
ImageFlattenProcessOutputDir = funParams.OutputDirectory;

if (~exist(ImageFlattenProcessOutputDir,'dir'))
    mkdir(ImageFlattenProcessOutputDir);
end

nFrame = movieData.nFrames_;

for iChannel = selected_channels
    
    img_pixel_pool = [];
    for iFrame = 1 : nFrame
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        img_pixel_pool = [img_pixel_pool currentImg(:)];
    end
    
    img_min = min(img_pixel_pool(:));
    img_max = min(max(img_pixel_pool(:)),  mean(img_pixel_pool(:))+10*std(img_pixel_pool(:)));
    
    [hist_all_frame, hist_bin] = hist(img_pixel_pool,25);
    
    for iFrame = 1 : nFrame
        hist_this_frame = hist_all_frame(:,iFrame);
        ind = find(hist_this_frame==max(hist_this_frame));
        center_value(iFrame) = hist_bin(ind(1));
    end
    
    center_value = center_value/max(center_value);
    center_value = sqrt(center_value);
    center_value = imfilter(center_value,[1 2 3 9 3 2 1]/21,'replicate','same');
    center_value = center_value/max(center_value);
    
    % Make output directory for the flattened images
    ImageFlattenChannelOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(ImageFlattenChannelOutputDir,'dir'))
        mkdir(ImageFlattenChannelOutputDir);
    end
    
    movieData.processes_{indexFlattenProcess}.setOutImagePath(iChannel,ImageFlattenChannelOutputDir)

    display(['Start to do image flatten in Channel ',num2str(iChannel)]);

    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImg = double(currentImg);
        
        % based on the given method index, do log or sqrt to flatten the image
        if flatten_method_ind == 1
            currentImg = log(currentImg);
            currentImg = currentImg - log(img_min);
            currentImg = currentImg/log((center_value(iFrame))*img_max);
            
        else
            if flatten_method_ind == 2
                
                currentImg = currentImg - img_min;
                currentImg = currentImg/(center_value(iFrame))/img_max;
    
                currentImg = sqrt(currentImg);
            end
        end

        % Smooth the image in requested
        if Gaussian_sigma > 0
            currentImg = imfilter(currentImg, fspecial('gaussian',round(5*Gaussian_sigma), Gaussian_sigma),'replicate','same');
        end
        currentImg(find(currentImg>1)) = 1;
        imwrite(currentImg,[ImageFlattenChannelOutputDir,'/flatten_',num2str(iFrame),'.tif']);
    end
end
