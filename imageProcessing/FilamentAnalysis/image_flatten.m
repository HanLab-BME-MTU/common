function movieData = image_flatten(movieData, varargin)

% Find the package of Filament Analysis
nPackage = length(movieData.packages_);

indexFilamentPackage = 0;
for i = 1 : nPackage
    if(strcmp(movieData.packages_{i}.getName,'FilamentAnalysis')==1)
        indexFilamentPackage = i;
        break;
    end
end

if(indexFilamentPackage==0)
    msg('Need to be in Filament Package for now.')
    return; 
end


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

TimeFilterSigma = funParams.TimeFilterSigma;                        
Sub_Sample_Num  = funParams.Sub_Sample_Num;

nFrame = movieData.nFrames_;


ImageFlattenProcessOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'ImageFlatten'];
if (~exist(ImageFlattenProcessOutputDir,'dir'))
    mkdir(ImageFlattenProcessOutputDir);
end

for iChannel = selected_channels
    ImageFlattenChannelOutputDir = [ImageFlattenProcessOutputDir,'/Channel',num2str(iChannel)];
    if (~exist(ImageFlattenChannelOutputDir,'dir'))
        mkdir(ImageFlattenChannelOutputDir);
    end
    
    movieData.processes_{indexFlattenProcess}.setOutImagePath(iChannel,ImageFlattenChannelOutputDir);
end

for iChannel = selected_channels

    ImageFlattenProcessOutputDir = movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel};
    
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    filename_short_strs = uncommon_str_takeout(movieData.channels_(iChannel).fileNames_);
    
    Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
    Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
    Frames_results_correspondence = Frames_results_correspondence(1:nFrame);
       
    img_pixel_pool = [];
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        img_pixel_pool = [img_pixel_pool currentImg(:)];
    end
    
    img_min = min(img_pixel_pool(:));
    img_max = min(max(img_pixel_pool(:)),  mean(img_pixel_pool(:))+10*std(img_pixel_pool(:)));
    
    [hist_all_frame, hist_bin] = hist(img_pixel_pool,55);
    
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        hist_this_frame = hist_all_frame(:,iFrame_subsample);
        ind = find(hist_this_frame==max(hist_this_frame));
        center_value(iFrame_subsample) = hist_bin(ind(1));
    end
    
    center_value = center_value/max(center_value);
    center_value = sqrt(center_value);
    center_value = imfilter(center_value,[1 2 3 9 3 2 1]/21,'replicate','same');
    center_value = center_value/max(center_value);
    
    % Make output directory for the flattened images
    ImageFlattenChannelOutputDir = movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel};
   if (~exist(ImageFlattenChannelOutputDir,'dir'))
        mkdir(ImageFlattenChannelOutputDir);
    end
   
    display(['Start to do image flatten in Channel ',num2str(iChannel)]);

     for iFrame_subsample = 1 : length(Frames_to_Seg)
         iFrame = Frames_to_Seg(iFrame_subsample);
         disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImg = double(currentImg);
        
        % based on the given method index, do log or sqrt to flatten the image
        if flatten_method_ind == 1
            currentImg = log(currentImg);
            currentImg = currentImg - log(img_min);
            currentImg = currentImg/log((center_value(iFrame_subsample))*img_max);
            
        else
            
            currentImg = currentImg - img_min;
            currentImg = currentImg/(center_value(iFrame_subsample))/img_max;
            if flatten_method_ind == 2
                currentImg = (currentImg).^(1/2);
            end
            
            if flatten_method_ind == 3
                currentImg = (currentImg).^(2/3);
            end
        end
    

        % Smooth the image in requested
        if Gaussian_sigma > 0
            currentImg = imfilter(currentImg, fspecial('gaussian',round(5*Gaussian_sigma), Gaussian_sigma),'replicate','same');
        end
        currentImg(find(currentImg>1)) = 1;
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame

              imwrite(currentImg,[ImageFlattenChannelOutputDir,'/flatten_', ...
                filename_short_strs{iFrame + sub_i-1},'.tif']);
           
            end
        end
        
        
        if(TimeFilterSigma > 0)
            if iFrame_subsample==1
                Image_tensor = zeros(size(currentImg,1),size(currentImg,2),length(Frames_to_Seg));
            end
            Image_tensor(:,:,iFrame_subsample) = currentImg;
        end
    end
end

if(TimeFilterSigma > 0)    
    FilterHalfLength = 2*ceil(TimeFilterSigma);
    
    temperal_filter = zeros(1,1,2*FilterHalfLength+1);
    
    H = fspecial('gaussian',2*FilterHalfLength+1, TimeFilterSigma);
    H_1D = H(FilterHalfLength+1,:);
    
    H_1D = H_1D/(sum(H_1D));
    
    temperal_filter(1,1,:) = H_1D(:);
    
    time_filtered = imfilter(Image_tensor,temperal_filter,'replicate','same');
    
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        disp(['Frame: ',num2str(iFrame)]);
        currentImg = squeeze(time_filtered(:,:,iFrame_subsample));
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
            disp(['Frame: ',num2str(iFrame + sub_i-1)]);
     
            imwrite(currentImg,[ImageFlattenChannelOutputDir,'/flatten_', ...
                filename_short_strs{iFrame + sub_i-1},'.tif']);
            end
        end
    end
end