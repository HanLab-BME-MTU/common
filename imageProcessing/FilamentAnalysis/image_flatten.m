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
    
    % Make output directory for the flattened images
    ImageFlattenChannelOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(ImageFlattenChannelOutputDir,'dir'))
        mkdir(ImageFlattenChannelOutputDir);
    end
    
    movieData.processes_{indexFlattenProcess}.setOutImagePath(iChannel,ImageFlattenChannelOutputDir)
    
    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImg = double(currentImg);
        
        % Normalize before flatten
%         currentImg = currentImg - min(min(currentImg)) +1;
        currentImg = currentImg/(max(max(currentImg)));
        
        % based on the given method index, do log or sqrt to flatten the image
        if flatten_method_ind == 1
            currentImg = log(currentImg);
            currentImg = currentImg - min(min(currentImg));
            currentImg = currentImg/(max(max(currentImg)));
            
        else
            if flatten_method_ind == 2
                currentImg = sqrt(currentImg);
            end
        end

        % Smooth the image in requested
        if Gaussian_sigma > 0
            currentImg = imfilter(currentImg, fspecial('gaussian',round(5*Gaussian_sigma), Gaussian_sigma),'replicate','same');
        end
        imwrite(currentImg,[ImageFlattenChannelOutputDir,'/flatten_',num2str(iFrame),'.tif']);
    end
end
