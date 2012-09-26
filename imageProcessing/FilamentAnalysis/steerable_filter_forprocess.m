function steerable_filter_forprocess(movieData, varargin)
% Created 07 2012 by Liya Ding, Matlab R2011b

nProcesses = length(movieData.processes_);

indexSteerabeleProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Steerable filtering')==1)
        indexSteerabeleProcess = i;
        break;
    end
end

if indexSteerabeleProcess==0
    msg('Please set parameters for steerable filtering.')
    return;
end
funParams=movieData.processes_{indexSteerabeleProcess}.funParams_;

selected_channels = funParams.ChannelIndex;
BaseSteerableFilterSigma = funParams.BaseSteerableFilterSigma;
Levelsofsteerablefilters = funParams.Levelsofsteerablefilters;
ImageFlattenFlag = funParams.ImageFlattenFlag;
Sub_Sample_Num = funParams.Sub_Sample_Num;

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess==0 && ImageFlattenFlag == 2
    display('The setting shows you want to use flattened image for steerable filtering. Please set parameters for Image Flatten and run.')
    return;
end


ImageSteerableFilterProcessOutputDir = funParams.OutputDirectory;
if (~exist(ImageSteerableFilterProcessOutputDir,'dir'))
    mkdir(ImageSteerableFilterProcessOutputDir);
end

nFrame = movieData.nFrames_;

% RES_cell = cell(1,nFrame);
% Image_cell = cell(1,nFrame);

Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
Frames_results_correspondence = Frames_results_correspondence(1:nFrame);

for iChannel = selected_channels
        % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    filename_short_strs = uncommon_str_takeout(movieData.channels_(iChannel).fileNames_);
    
    % Make output directory for the steerable filtered images
    ImageSteerableFilterChannelOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(ImageSteerableFilterChannelOutputDir,'dir'))
        mkdir(ImageSteerableFilterChannelOutputDir);
    end
    
%     if indexFlattenProcess >0
%         FileNames = movieData.processes_{indexFlattenProcess}.getOutImageFileNames(iChannel);
%     end
%     
    display(['Start to do steerable filtering in Channel ',num2str(iChannel)]);

   
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        currentImg = double(currentImg);

        levels_sizes = 2.^((1:Levelsofsteerablefilters)-1);
        
        % Steerable filtering using four scales one doubling the previous one.
        % function multiscaleSteerableDetector will automatically merge the results
        [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
       
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite((MAX_st_res)/(max(max(MAX_st_res))), ...
                    [ImageSteerableFilterChannelOutputDir,'/MAX_st_res_', ...
                    filename_short_strs{iFrame + sub_i-1},'.tif']);
            end
        end
        
        save([ImageSteerableFilterChannelOutputDir,'/steerable_',filename_short_strs{iFrame},'.mat'],...
            'orienation_map', 'MAX_st_res','nms');
    end
end
