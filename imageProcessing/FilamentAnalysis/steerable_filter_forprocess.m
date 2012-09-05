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
% 
%         if iFrame==1
%             Image_tensor = zeros(size(currentImg,1),size(currentImg,2),nFrame);
%             Res_tensor = zeros(size(currentImg,1),size(currentImg,2),nFrame);
%         end
        currentImg = double(currentImg);
%         Image_cell{iFrame} = currentImg;
%         Image_tensor(:,:,iFrame) = currentImg;
%         levels = 0:Levelsofsteerablefilters-1;
        levels_sizes = 2.^((1:Levelsofsteerablefilters)-1);
        
        % Steerable filtering using four scales one doubling the previous one.
        % function multiscaleSteerableDetector will automatically merge the results
         [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
%        [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma+levels);
        
        imwrite((MAX_st_res)/(max(max(MAX_st_res))),[ImageSteerableFilterChannelOutputDir,'/MAX_st_res_',num2str(iFrame),'.tif']);
        
        save([ImageSteerableFilterChannelOutputDir,'/steerable_',num2str(iFrame),'.mat'],...
            'orienation_map', 'MAX_st_res','nms');
        
%         RES_cell{iFrame} = MAX_st_res;
%          Res_tensor(:,:,iFrame) = MAX_st_res;
    end
    
end
% 
% save('temperal.mat','Image_cell','RES_cell','Image_tensor','Res_tensor');
% 
% temperal_filter = zeros(1,1,11);
% 
% H = fspecial('gaussian',11, 2);
% H_1D = H(6,:);
% 
% H_1D = H_1D/(sum(H_1D));
% 
% temperal_filter(1,1,:) = H_1D(:);
% 
% time_filtered = imfilter(Image_tensor,temperal_filter,'replicate','same');
%         
% 
% 
%     for iFrame = 1 : nFrame
%         disp(['Frame: ',num2str(iFrame)]);
%         
%         currentImg = squeeze(time_filtered(:,:,iFrame));
%         
%         levels_sizes = 2.^((1:Levelsofsteerablefilters)-1);
%         
%         % Steerable filtering using four scales one doubling the previous one.
%         % function multiscaleSteerableDetector will automatically merge the results
%          [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
% %        [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma+levels);
%         
%         imwrite([RES_cell{iFrame}/(max(max(RES_cell{iFrame}))) (MAX_st_res)/(max(max(RES_cell{iFrame})))],[ImageSteerableFilterChannelOutputDir,'/time_filtered_MAX_st_res_',num2str(iFrame),'.tif']);
%         display([mean2(RES_cell{iFrame})   mean2(MAX_st_res)]);
%         save([ImageSteerableFilterChannelOutputDir,'/time_filtered_steerable_',num2str(iFrame),'.mat'],...
%             'orienation_map', 'MAX_st_res','nms');        
%     end