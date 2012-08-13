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

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess==0
    display('If you need flattening, please set parameters for Image Flatten and run. For now, we use the original images.')
%     return;
end

funParams=movieData.processes_{indexSteerabeleProcess}.funParams_;

selected_channels = funParams.ChannelIndex;
BaseSteerableFilterSigma = funParams.BaseSteerableFilterSigma;
Levelsofsteerablefilters = funParams.Levelsofsteerablefilters;

ImageSteerableFilterProcessOutputDir = funParams.OutputDirectory;
if (~exist(ImageSteerableFilterProcessOutputDir,'dir'))
    mkdir(ImageSteerableFilterProcessOutputDir);
end

nFrame = movieData.nFrames_;

for iChannel = selected_channels
    
    % Make output directory for the steerable filtered images
    ImageSteerableFilterChannelOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(ImageSteerableFilterChannelOutputDir,'dir'))
        mkdir(ImageSteerableFilterChannelOutputDir);
    end
    
    if indexFlattenProcess >0
        FileNames = movieData.processes_{indexFlattenProcess}.getOutImageFileNames(iChannel);
    end
    
    display(['Start to do steerable filtering in Channel ',num2str(iChannel)]);

    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, FileNames{1}{iFrame}]);
        else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        
        currentImg = double(currentImg);
        
%         levels = 0:Levelsofsteerablefilters-1;
        levels_sizes = 2.^(Levelsofsteerablefilters-1);
        
        % Steerable filtering using four scales one doubling the previous one.
        % function multiscaleSteerableDetector will automatically merge the results
         [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
%        [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma+levels);
        
        imwrite((MAX_st_res)/(max(max(MAX_st_res))),[ImageSteerableFilterChannelOutputDir,'/MAX_st_res_',num2str(iFrame),'.tif']);
        
        save([ImageSteerableFilterChannelOutputDir,'/steerable_',num2str(iFrame),'.mat'],...
            'orienation_map', 'MAX_st_res');
    end
end
