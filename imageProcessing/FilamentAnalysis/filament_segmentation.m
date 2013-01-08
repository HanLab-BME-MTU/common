function movieData = filament_segmentation(movieData, varargin)
% Created 07 2012 by Liya Ding, Matlab R2011b

% input movieData object, with the parameters
%   funParams.ChannelIndex:         the channels to process
%   funParams.Pace_Size:            the parameter to set pace in local segmentation
%   funParams.Patch_Size:           the parameter to set patch size in local segmentation, for the estimation of local threshold
%   funParams.lowerbound_localthresholding: The percentage as the lower bound of local thresholding
%                                    local threshold has to be larger or equal to this percentage of the global threshold
%   funParams.Combine_Way :         The way to combine segmentation results from steerable filtering responce
%                                     and from intensity, default is : only use steerable filtering result
%   funParams.Cell_Mask_ind:        Flag to set if cell mask is used, if 1, use segmentation(refined) results,
%                                     if 2, use the user define ROI as in MD_ROI.tif in movieData folder, if 3, no such limit
%   funParams.VIF_Outgrowth_Flag:   Flag to do VIF_outgrowth or not. This is an option made for Gelfand lab

% a temp flag for saving tif stack image for Gelfand Lab
save_tif_flag=1;


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


nProcesses = length(movieData.processes_);

indexFilamentSegmentationProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Filament Segmentation')==1)
        indexFilamentSegmentationProcess = i;
        break;
    end
end

if indexFilamentSegmentationProcess==0
    msg('Please set parameters for steerable filtering.')
    return;
end


funParams=movieData.processes_{indexFilamentSegmentationProcess}.funParams_;

selected_channels = funParams.ChannelIndex;
StPace_Size = funParams.StPace_Size;
StPatch_Size = funParams.StPatch_Size;
Stlowerbound =  funParams.st_lowerbound_localthresholding;
IntPace_Size = funParams.IntPace_Size;
IntPatch_Size = funParams.IntPatch_Size;
Intlowerbound =  funParams.int_lowerbound_localthresholding;

Combine_Way = funParams.Combine_Way;
Cell_Mask_ind = funParams.Cell_Mask_ind;
VIF_Outgrowth_Flag = funParams.VIF_Outgrowth_Flag;
Sub_Sample_Num  = funParams.Sub_Sample_Num;

%% Output Directories

FilamentSegmentationProcessOutputDir  = [movieData.packages_{indexFilamentPackage}.outputDirectory_, filesep 'FilamentSegmentation'];
if (~exist(FilamentSegmentationProcessOutputDir,'dir'))
    mkdir(FilamentSegmentationProcessOutputDir);
end

for iChannel = selected_channels
    FilamentSegmentationChannelOutputDir = [FilamentSegmentationProcessOutputDir,'/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    movieData.processes_{indexFilamentSegmentationProcess}.setOutImagePath(iChannel,FilamentSegmentationChannelOutputDir);
end



%%
indexSteerabeleProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Steerable filtering')==1)
        indexSteerabeleProcess = i;
        break;
    end
end

if indexSteerabeleProcess==0 && Combine_Way~=2
    msg('Please run steerable filtering first.')
    return;
end

funParams_st=movieData.processes_{indexSteerabeleProcess}.funParams_;

BaseSteerableFilterSigma = funParams_st.BaseSteerableFilterSigma;
Levelsofsteerablefilters = funParams_st.Levelsofsteerablefilters;
ImageFlattenFlag = funParams_st.ImageFlattenFlag;

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

if indexFlattenProcess == 0 && ImageFlattenFlag==2
    display('Please set parameters for Image Flatten.')
    return;
end

indexCellSegProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Mask Refinement')==1)
        indexCellSegProcess = i;
        break;
    end
end

if indexCellSegProcess == 0 && Cell_Mask_ind == 1
    msg('Please run segmentation and refinement first.')
    return;
end


nFrame = movieData.nFrames_;

% If the user set an cell ROI read in
if(exist([movieData.outputDirectory_,filesep,'MD_ROI.tif'],'file'))
    user_input_mask = imread([movieData.outputDirectory_,filesep,'MD_ROI.tif']);
end

%% Prepare the cone masks
cone_size = 15;
cone_angle = 25;
cone_mask = cell(180,1);
cone_zero = zeros(2*cone_size+1,2*cone_size+1);
for ci = 1 : 2*cone_size+1
    for cj = 1 : 2*cone_size+1
        cone_zero(ci,cj) = mod(atan2(ci-cone_size-1,cj-cone_size-1),pi);
    end
end
cone_zero_mask = cone_zero<cone_angle/180*pi | cone_zero>pi - cone_angle/180*pi;

cone_zero_mask(cone_size-3:cone_size+3,:)=1;


for cone_i = 1 :180
    cone_mask{cone_i} = imrotate(cone_zero_mask, cone_i, 'nearest','crop');
end

%%
for iChannel = selected_channels
    
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
   Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir =  movieData.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel};
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    HeatOutputDir = [FilamentSegmentationChannelOutputDir,'/HeatOutput'];
    
    if (~exist(HeatOutputDir,'dir'))
        mkdir(HeatOutputDir);
    end
    
    HeatEnhOutputDir = [HeatOutputDir,'/Enh'];
    
    if (~exist(HeatEnhOutputDir,'dir'))
        mkdir(HeatEnhOutputDir);
    end
    
    DataOutputDir = [FilamentSegmentationChannelOutputDir,'/DataOutput'];
    
    if (~exist(DataOutputDir,'dir'))
        mkdir(DataOutputDir);
    end
    
    
    
    
    % If steerable filter process is run
    if indexSteerabeleProcess>0
        SteerableChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    end
    
    %     if indexFlattenProcess >0
    %         FileNames = movieData.processes_{indexFlattenProcess}.getOutImageFileNames(iChannel);
    %     end
    %
    display('======================================');
    display(['Current movie: as in ',movieData.outputDirectory_]);
    display(['Start filament segmentation in Channel ',num2str(iChannel)]);
    
    % Segment only the real collected data, but skip the padded ones, which
    % were there just to fill in the time lap to make two channel same
    % number of frames
    Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
    Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
    Frames_results_correspondence = Frames_results_correspondence(1:nFrame);
    
%     indexFlattenProcess=1;
    for iFrame_index = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_index);
        
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0 && ImageFlattenFlag==2
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        currentImg = double(currentImg);
        
        
        if( save_tif_flag==1 && iFrame==Frames_to_Seg(1)  )
            % Gelfand lab needs single file results for tif stack
            tif_stack_binary_seg_image_data = uint8(zeros(size(currentImg,1),size(currentImg,2),length(Frames_to_Seg)));
            tif_stack_RGB_heat_image_data = uint8(zeros(size(currentImg,1),size(currentImg,2),3,length(Frames_to_Seg)));
        end
        
        load([SteerableChannelOutputDir, filesep, 'steerable_', ...
            filename_short_strs{iFrame},'.mat']);
        
        NMS_Segment=[];
        Intensity_Segment=[];
        SteerabelRes_Segment=[];
        
        Min_area = 20;
        Min_longaxis = 6;
        
        
        switch Combine_Way
            case 'int_st_both'
                level0 = thresholdOtsu(MAX_st_res);
                thresh_Segment = MAX_st_res > level0;
                
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,0);
                [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,0);
                current_seg = and(Intensity_Segment,SteerabelRes_Segment);
                
            case 'st_only'
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,0);
                current_seg = SteerabelRes_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                
            case 'st_nms_two'
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,0);
                [level2, NMS_Segment ] = thresholdLocalSeg(nms,'Rosin',StPatch_Size,StPace_Size,Stlowerbound,0);
                current_seg = imdilateWithScale(NMS_Segment,scaleMap,BaseSteerableFilterSigma.*(2.^((1:Levelsofsteerablefilters)-1)))...
                    .*SteerabelRes_Segment;
                       
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
        
            case 'st_nms_only'
                [level2, NMS_Segment ] = thresholdLocalSeg(nms,'Rosin',StPatch_Size,StPace_Size,Stlowerbound,0);
                current_seg = NMS_Segment;
                Intensity_Segment = current_seg;
                SteerabelRes_Segment = current_seg;
                Min_area = 6;
                
            case 'int_only'
                [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,0);
                
                current_seg = Intensity_Segment;
                SteerabelRes_Segment = current_seg;
            otherwise
                warning('Use the default of union');
                [level1, SteerabelRes_Segment ] = thresholdLocalSeg(MAX_st_res,'Otsu',StPatch_Size,StPace_Size,Stlowerbound,0);
                [level2, Intensity_Segment ] = thresholdLocalSeg(currentImg,'Otsu',IntPatch_Size,IntPace_Size,Intlowerbound,0);
                % The segmentation is set as the union of two segmentation.
                current_seg = or(Intensity_Segment,SteerabelRes_Segment);
        end
        
        if Cell_Mask_ind == 1 % using cell segmentation from same channel
            MaskCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame);
        else
            if Cell_Mask_ind == 2 % Using input static ROI tiff
                MaskCell = user_input_mask>0;
            else
                if Cell_Mask_ind == 4 % No limit
                    MaskCell = ones(size(currentImg,1),size(currentImg,2));
                else
                    % Combine from both channel
                    % In this option, the channel need to be 1. MT or Membrame, 2. VIF or Actin
                    MaskVIFCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(2,iFrame);
                    MaskMTCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(1,iFrame);
                    
                    H_close_cell = fspecial('disk',5);
                    H_close_cell = H_close_cell>0;
                    
                    MaskMTCell = imerode(MaskMTCell,H_close_cell);
                    TightMask = MaskMTCell.*MaskMTCell;
                    
                    % Make the mask bigger in order to include all
                    MaskCell = imdilate(TightMask, ones(15,15),'same');
                end
            end
        end
        
        current_seg = current_seg.*MaskCell;
        
        
        %%
        % A smoothing done only at the steerable filtering results, if only intensity only, then the same
        orienation_map_filtered = OrientationSmooth(orienation_map, SteerabelRes_Segment);
        
        %%
        % Voting of the orientation field for the non-steerable filter
        % segmented places.
        
        OrientationVoted = OrientationVote(orienation_map,SteerabelRes_Segment,3,45);
        
        intensity_addon = current_seg - SteerabelRes_Segment ==1;
        if (~isempty(max(max(intensity_addon))>0))
            orienation_map_filtered(find(intensity_addon>0)) = OrientationVoted(find(intensity_addon>0));
        end
        
        
                %% Deleting the small isolated dots
        
                labelMask = bwlabel(current_seg);
        
                ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','MinorAxisLength');
        
                obAreas = [ob_prop.Area];
                obLongaxis = [ob_prop.MajorAxisLength];
                obShortaxis = [ob_prop.MinorAxisLength];
                obEccentricity = [ob_prop.Eccentricity];
                ratio  = obShortaxis./obLongaxis;
                % for now the parameters are set here, without adaptiveness
                for i_area = 1 : length(obAreas)
                    if obAreas(i_area) < 100
                        angle_area{i_area} = orienation_map_filtered(find(labelMask==i_area));
                        [h_area, bin] = hist(angle_area{i_area},-pi/2:5/180*pi:pi/2);
                        ind_t = find(h_area==max(h_area));
                        temp = mod((angle_area{i_area} - bin(ind_t(1)) + pi/2), pi) - pi/2;
                        if std(temp)>0.75 && max(h_area)<0.2*length(angle_area{i_area}) && ratio(i_area) >0.5 && obLongaxis(i_area)<20
                            labelMask(find(labelMask==i_area))=0;
                        end
                    end
                end
        
                current_seg = labelMask > 0;
                
                
                       labelMask = bwlabel(current_seg);
        
        ob_prop = regionprops(labelMask,'Area','MajorAxisLength','Eccentricity','Centroid');
        
        if length(ob_prop) > 1
            obAreas = [ob_prop.Area];
            obLongaxis = [ob_prop.MajorAxisLength];
            obEccentricity = [ob_prop.Eccentricity];
            obCentroid = [ob_prop.Centroid];
            
            for i_area = 1 : length(obAreas)
                centroid_x = round(obCentroid(2*i_area-1));
                centroid_y = round(obCentroid(2*i_area));
                
                    if obAreas(i_area) <Min_area || obLongaxis(i_area) <Min_longaxis
                        labelMask(labelMask==i_area) = 0;
                    end
               
            end
        end
        
        current_seg = labelMask > 0;
        
        
        %
        %         [ind_a,ind_b] = find(current_seg>0);
        %
        %         cone_bins = cell(size(current_seg,1), size(current_seg,2));
        %         cone_weight_bins = cell(size(current_seg,1), size(current_seg,2));
        %
        %         for si = 1 : length(ind_a)
        %             pixel_angle = round((-orienation_map(ind_a(si), ind_b(si))+pi/2)*180/pi);
        %             weight_res = MAX_st_res(ind_a(si), ind_b(si));
        %             if pixel_angle ==0
        %                 pixel_angle = 180;
        %             end
        %             try
        %                 [ind_c,ind_d] = find(cone_mask{pixel_angle}>0);
        %                 for p_i = 1 : length(ind_c)
        %                     cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} ...
        %                         = [cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} pixel_angle];
        %                     cone_weight_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} ...
        %                         = [cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} weight_res];
        %                 end
        %             end
        %
        %         end
        %
        %         for p_i = 1 : size(cone_bins,1)
        %             for p_j = 1 : size(cone_bins,2)
        %                 length_cone_bin(p_i,p_j) = length(cone_bins{p_i,p_j});
        %                  weight_bin(p_i,p_j) = m(cone_weight_bins{p_i,p_j});
        %                 if(~isempty(cone_bins{p_i,p_j}))
        %                 h = hist(cone_bins{p_i,p_j},0:10:180);
        %                 centernumber_cone_bin(p_i,p_j) = max(h);
        %                 else
        %                      centernumber_cone_bin(p_i,p_j) =0;
        %                 end
        %             end
        %         end
        
        
        %% For heat presentation of the segmented filaments
        
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite(current_seg, ...
                    [FilamentSegmentationChannelOutputDir,'/segment_binary_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);
            end
        end
        
        
        currentImg = uint8(currentImg/1);
        Hue = (-orienation_map_filtered(:)+pi/2)/(pi)-0.2;
        Hue(find(Hue>=1)) = Hue(find(Hue>=1)) -1;
        Hue(find(Hue<0)) = Hue(find(Hue<0)) +1;
        
        Sat = Hue*0+1;
        Value = Hue*0+1;
        RGB_seg_orient_heat_array = hsv2rgb([Hue Sat Value]);
        R_seg_orient_heat_map = col2im(RGB_seg_orient_heat_array(:,1),[1 1],[size(current_seg,1) size(current_seg,2)]);
        G_seg_orient_heat_map = col2im(RGB_seg_orient_heat_array(:,2),[1 1],[size(current_seg,1) size(current_seg,2)]);
        B_seg_orient_heat_map = col2im(RGB_seg_orient_heat_array(:,3),[1 1],[size(current_seg,1) size(current_seg,2)]);
        
        enhanced_im_r = currentImg;
        enhanced_im_g = currentImg;
        enhanced_im_b = currentImg;
        
        enhanced_im_r(find(current_seg>0))=255*R_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_g(find(current_seg>0))=255*G_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_b(find(current_seg>0))=255*B_seg_orient_heat_map(find(current_seg>0));
        
        RGB_seg_orient_heat_map(:,:,1 ) = enhanced_im_r;
        RGB_seg_orient_heat_map(:,:,2 ) = enhanced_im_g;
        RGB_seg_orient_heat_map(:,:,3 ) = enhanced_im_b;
        
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite(RGB_seg_orient_heat_map, ...
                    [HeatEnhOutputDir,'/segment_heat_',...
                    filename_short_strs{iFrame+ sub_i-1},'.tif']);
            end
        end
        
        %% Save segmentation results
        save([DataOutputDir,'/steerable_vote_', ...
            filename_short_strs{iFrame},'.mat'],...
            'currentImg','orienation_map_filtered','OrientationVoted','orienation_map', ...
            'MAX_st_res', 'current_seg','Intensity_Segment','SteerabelRes_Segment','NMS_Segment');
        
        if( save_tif_flag==1)
%             current_seg = (imread([FilamentSegmentationChannelOutputDir,'/segment_binary_',filename_short_strs{iFrame},'.tif']))>0;
%             RGB_seg_orient_heat_map = imread([HeatEnhOutputDir,'/segment_heat_',filename_short_strs{iFrame},'.tif']);
%             
            tif_stack_binary_seg_image_data(:,:,iFrame_index) = uint8(current_seg*255);
            tif_stack_RGB_heat_image_data(:,:,:,iFrame_index) = uint8(RGB_seg_orient_heat_map);
            
        end
    end
    %% For Gelfand Lab, save results as tif stack file
    if( save_tif_flag==1)
        options.comp = false;
        options.ask = false;
        options.message = true;
        options.append = false;
        
        % Save the multi-frame RGB color image
        options.color = true;
        saveastiff(tif_stack_RGB_heat_image_data, [FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_seg_heat.tif'], options);
        options.color = false;
        saveastiff(tif_stack_binary_seg_image_data, [FilamentSegmentationProcessOutputDir,'/channel_',num2str(iChannel),'_seg_binary.tif'], options);
        
    end
    
end



%% For Gelfand Lab, outgrowth calculation
if(VIF_Outgrowth_Flag==1)
    VIF_outgrowth_measurement(movieData);
end




