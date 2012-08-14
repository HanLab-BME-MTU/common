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


indexSteerabeleProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Steerable filtering')==1)
        indexSteerabeleProcess = i;
        break;
    end
end

if indexSteerabeleProcess==0
    msg('Please run steerable filtering first.')
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
    display('Please set parameters for Image Flatten.')
%     return;
end

indexCellSegProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Mask Refinement')==1)
        indexCellSegProcess = i;
        break;
    end
end

if indexCellSegProcess==0
    msg('Please run segmentation and refinement first.')
    return;
end

funParams=movieData.processes_{indexFilamentSegmentationProcess}.funParams_;

selected_channels = funParams.ChannelIndex;
Pace_Size = funParams.Pace_Size;
Patch_Size = funParams.Patch_Size;
Combine_Way = funParams.Combine_Way;
Cell_Mask_ind = funParams.Cell_Mask_ind;
lowerbound =  funParams.lowerbound_localthresholding;
VIF_Outgrowth_Flag = funParams.VIF_Outgrowth_Flag;

FilamentSegmentationOutputDir = funParams.OutputDirectory;

if (~exist(FilamentSegmentationOutputDir,'dir'))
    mkdir(FilamentSegmentationOutputDir);
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
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end

    SteerableChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};

    
    if indexFlattenProcess >0
        FileNames = movieData.processes_{indexFlattenProcess}.getOutImageFileNames(iChannel);
    end
    
    display(['Start to do filament segmentation in Channel ',num2str(iChannel)]);

    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, FileNames{1}{iFrame}]);
         else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        currentImg = double(currentImg);
        
        load([SteerableChannelOutputDir, filesep, 'steerable_',num2str(iFrame),'.mat']);
        
        MaskCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame);
       
        switch Combine_Way
            case 'int_st_both'
                 level0 = thresholdOtsu(MAX_st_res);
                 thresh_Segment = MAX_st_res > level0;
                                
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,Patch_Size,Pace_Size,lowerbound,0);
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,Patch_Size,Pace_Size,lowerbound,0);
                
            case 'st_only'                
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,Patch_Size,Pace_Size,lowerbound,0);
                current_seg = SteerabelRes_Segment; 
                Intensity_Segment = current_seg*0;
                
            case 'int_only'
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,Patch_Size,Pace_Size,lowerbound,0);
                current_seg = Intensity_Segment; 
                SteerabelRes_Segment = current_seg*0;
            otherwise
                warning('Use the default of union');
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,Patch_Size,Pace_Size,lowerbound,0);
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,Patch_Size,Pace_Size,lowerbound,0);
                % The segmentation is set as the union of two segmentation.
                current_seg = or(Intensity_Segment,SteerabelRes_Segment);
        end
        
        if Cell_Mask_ind == 1
            MaskCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame);
        else
            if Cell_Mask_ind == 2
                MaskCell = user_input_mask>0;
            else
                MaskCell = ones(size(currentImg,1),size(currentImg,2));
            end
        end
        
        current_seg = current_seg.*MaskCell;        

        
        %%
        % A smoothing done only at the steerable filtering results        
        orienation_map_filtered = OrientationSmooth(orienation_map, SteerabelRes_Segment);
        
        %%
        % Voting of the orientation field for the non-steerable filter
        % segmented places.
        tic
        OrientationVoted = OrientationVote(orienation_map,SteerabelRes_Segment,3,25);
        toc
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
            if obAreas(i_area) < 200
                angle_area{i_area} = orienation_map_filtered(find(labelMask==i_area));
                [h_area, bin] = hist(angle_area{i_area},-pi/2:5/180*pi:pi/2);
                ind_t = find(h_area==max(h_area));
                temp = mod((angle_area{i_area} - bin(ind_t(1)) + pi/2), pi) - pi/2;
                if std(temp>0.75) && max(h_area)<0.2*length(angle_area{i_area}) && ratio(i_area) >0.5 && obLongaxis(i_area)<20
                    labelMask(find(labelMask==i_area))=0;
                end
            end
        end
        
        current_seg = labelMask > 0;
%     
%         [ind_a,ind_b] = find(current_seg>0);
%         
%         cone_bins = cell(size(current_seg,1), size(current_seg,2));
%         
%         for si = 1 : length(ind_a)
%             pixel_angle = round(orienation_map(ind_a(si), ind_b(si))*180/pi);
%             if pixel_angle ==0
%                 pixel_angle = 180;
%             end
%             try
%                 [ind_c,ind_d] = find(cone_mask{pixel_angle}>0);
%                 for p_i = 1 : length(ind_c)
%                     cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} ...
%                         = [cone_bins{ind_c(p_i)+ind_a(si)-cone_size-1, ind_d(p_i)+ind_b(si)-cone_size-1} pixel_angle];
%                 end
%             end
%             
%         end
%         
%         for p_i = 1 : size(cone_bins,1)
%             for p_j = 1 : size(cone_bins,2)
%                 length_cone_bin(p_i,p_j) = length(cone_bins{p_i,p_j});
%                 if(~isempty(cone_bins{p_i,p_j}))
%                 h = hist(cone_bins{p_i,p_j},0:10:180);
%                 centernumber_cone_bin(p_i,p_j) = max(h);
%                 end
%             end
%         end
%                    
        
        %% For heat presentation of the segmented filaments
        
%         imwrite(current_seg,[FilamentSegmentationChannelOutputDir,'/segment_',num2str(iFrame),'.tif']);
        currentImg = uint8(currentImg);
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
        
        imwrite(RGB_seg_orient_heat_map,[FilamentSegmentationChannelOutputDir,'/segment_heat_',num2str(iFrame),'.tif']);
        
        %% Save segmentation results
        save([FilamentSegmentationChannelOutputDir,'/steerable_vote_',num2str(iFrame),'.mat'],...
            'currentImg','orienation_map_filtered','OrientationVoted','orienation_map', ...
            'MAX_st_res', 'current_seg','Intensity_Segment','SteerabelRes_Segment');
        
    end
end     
 

%% For Gelfand Lab, outgrowth calculation
if(VIF_Outgrowth_Flag==1)
   VIF_outgrowth_measurement(movieData);
end


