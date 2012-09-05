function VIF_outgrowth_measurement(movieData)

% Created June 2012 by Liya Ding, Matlab R2011b

nProcesses = length(movieData.processes_);

indexFilamentSegmentationProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Filament Segmentation')==1)
        indexFilamentSegmentationProcess = i;
        break;
    end
end

if indexFilamentSegmentationProcess==0
    msg('Please set parameters for Filament Segmentation and run.')
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
else
    funParams_st = movieData.processes_{indexSteerabeleProcess}.funParams_;
    ImageFlattenFlag = funParams_st.ImageFlattenFlag;    
end

indexFlattenProcess = 0;
for i = 1 : nProcesses
    if(strcmp(movieData.processes_{i}.getName,'Image Flatten')==1)
        indexFlattenProcess = i;
        break;
    end
end

 if indexFlattenProcess==0  && ImageFlattenFlag == 2
    display('The setting shows you want to use flattened image for steerable filtering. Please set parameters for Image Flatten and run.')
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

funParams=movieData.processes_{indexFilamentSegmentationProcess}.funParams_;

FilamentSegmentationOutputDir = funParams.OutputDirectory;

if (~exist(FilamentSegmentationOutputDir,'dir'))
    mkdir(FilamentSegmentationOutputDir);
end

nFrame = movieData.nFrames_;
selected_channels = [1];

% If the user set an cell ROI read in
if(exist([movieData.outputDirectory_,filesep,'MD_ROI.tif'],'file'))
    user_input_mask = imread([movieData.outputDirectory_,filesep,'MD_ROI.tif']);
end


% Gets the screen size
scrsz = get(0,'ScreenSize');
save_fig_flag = 0;

for iChannel = selected_channels
    
    % Make output directory for the steerable filtered images
    FilamentSegmentationChannelOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationChannelOutputDir,'dir'))
        mkdir(FilamentSegmentationChannelOutputDir);
    end
    
    SteerableChannelOutputDir = movieData.processes_{indexSteerabeleProcess}.outFilePaths_{iChannel};
    
    HeatOutputDir = [FilamentSegmentationChannelOutputDir,'/HeatOutput'];
    
    if (~exist(HeatOutputDir,'dir'))
        mkdir(HeatOutputDir);
    end
    
    HeatEnhOutputDir = [HeatOutputDir,'/Enh'];
    
    if (~exist(HeatEnhOutputDir,'dir'))
        mkdir(HeatEnhOutputDir);
    end
    
    HeatEnhBoundOutputDir = [HeatOutputDir,'/Enh_bound'];
    
    if (~exist(HeatEnhBoundOutputDir,'dir'))
        mkdir(HeatEnhBoundOutputDir);
    end
    
    if indexFlattenProcess >0
        FileNames = movieData.processes_{indexFlattenProcess}.getOutImageFileNames(iChannel);
    end
    disp(funParams.OutputDirectory);
    
    H_close = fspecial('disk',39);
    H_close = H_close>0;
    

    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
         % Read in the intensity image, flattened or original
         if indexFlattenProcess > 0
             currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, FileNames{1}{iFrame}]);
         else
             currentImg = movieData.channels_(iChannel).loadImage(iFrame);
         end
        
        load([SteerableChannelOutputDir, filesep, 'steerable_',num2str(iFrame),'.mat']);
        
        if funParams.Cell_Mask_ind == 1
            MaskCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame);
        else
            if funParams.Cell_Mask_ind == 2
                MaskCell = user_input_mask>0;
            else
                MaskCell = ones(size(currentImg,1),size(currentImg,2));
            end
        end
        load([FilamentSegmentationChannelOutputDir,'/steerable_vote_',num2str(iFrame),'.mat'],...
            'orienation_map_filtered','OrientationVoted','orienation_map', ...
            'MAX_st_res', 'current_seg','Intensity_Segment','SteerabelRes_Segment');
        
        if iFrame==1
            % Read in the initial circle from the 'start_ROI.tif' file
            MaskFirstFrame = imread([movieData.outputDirectory_,'/start_ROI.tif']);
            MaskFirstFrame = (MaskFirstFrame)>0;
            
            RoiYX = bwboundaries(MaskFirstFrame);
            RoiYX = RoiYX{1};
        end
        
        current_seg = SteerabelRes_Segment.*MaskCell;
        %current_seg = SteerabelRes_Segment;
        
        incircle_seg = SteerabelRes_Segment.*MaskFirstFrame;
        
        incircle_seg_dilate = imdilate(incircle_seg, H_close);
         
        all_seg_dilate = imdilate(current_seg, H_close);
         
        ind_incircle = find(incircle_seg_dilate>0);
        labelMask = bwlabel(all_seg_dilate);
        dilate_keep = labelMask.*0;
        
        ind_new_labels = unique(labelMask(ind_incircle));
       
        for li = 1 :  length(ind_new_labels)
            dilate_keep(find(labelMask == ind_new_labels(li))) = ind_new_labels(li);
        end

%         labelMask(ind_incircle) = keep_largest_area(all_seg_dilate);

        current_seg = dilate_keep.*current_seg;

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
                
                if( MaskFirstFrame(centroid_y, centroid_x) == 0)
                    if obAreas(i_area) < 12 || obLongaxis(i_area) < 10
                        labelMask(labelMask==i_area) = 0;
                    end
                else
                    if obAreas(i_area) < 7 || obLongaxis(i_area) < 5
                        labelMask(labelMask==i_area) = 0;
                    end
                end
            end
        end
        
        current_seg = labelMask > 0;
        
        if iFrame==1
            % Get the first segmented results for the base of comparison
            current_seg_firstframe = current_seg;
        
            current_seg_inside_firstframe = current_seg_firstframe.*(double(MaskFirstFrame));
            seg_sum_inside_firstframe = sum(sum(current_seg_inside_firstframe));       
        end
        
        
        current_seg_outside = current_seg.*(1-MaskFirstFrame);
        
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
        
        enhanced_im_r(find(current_seg>0)) = 255*R_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_g(find(current_seg>0)) = 255*G_seg_orient_heat_map(find(current_seg>0));
        enhanced_im_b(find(current_seg>0)) = 255*B_seg_orient_heat_map(find(current_seg>0));
        
        RGB_seg_orient_heat_map(:,:,1) = enhanced_im_r;
        RGB_seg_orient_heat_map(:,:,2) = enhanced_im_g;
        RGB_seg_orient_heat_map(:,:,3) = enhanced_im_b;
        
%         RGB_seg_orient_heat_map = uint8(RGB_seg_orient_heat_map*255);
        
        h12 = figure(12);
        hold off;
        
        imagesc(RGB_seg_orient_heat_map);
        axis image; axis off;
        title(['Percentage of out growth: ', ...
            num2str(100*sum(sum(current_seg_outside))/seg_sum_inside_firstframe), '%'],'FontSize',20);
        
        saveas(h12,[HeatEnhOutputDir,'/Enh_VIF_heat_display_',num2str(iFrame),'.tif']);
        if(save_fig_flag==1)
            saveas(h12,[HeatEnhOutputDir,'/Enh_VIF_heat_display_',num2str(iFrame),'.fig']);
        end
        
        hold on; plot(RoiYX(:,2),RoiYX(:,1),'m');
        
        saveas(h12,[HeatEnhOutputDir,'_bound/Enh_Bound_VIF_heat_display_',num2str(iFrame),'.tif']);
        if(save_fig_flag==1)
            saveas(h12,[HeatEnhOutputDir,'_bound/Enh_Bound_VIF_heat_display_',num2str(iFrame),'.fig']);
        end
        seg_outside_current(iChannel, iFrame) = sum(sum(current_seg_outside));
        ratio_outside_firstframeinside(iChannel, iFrame) = sum(sum(current_seg_outside))/seg_sum_inside_firstframe;
        
    end
end

display('The outgrowth percentage results:');
display(ratio_outside_firstframeinside'*100);

% Save the outgrowth results
save([FilamentSegmentationChannelOutputDir,'/seg_outside.mat'],'seg_outside_current','seg_sum_inside_firstframe','ratio_outside_firstframeinside');




%% Simply use intensity as outgrowth calculation

for iFrame = 1 : nFrame
    disp(['Frame: ',num2str(iFrame)]);
    
    % Read in the intensity image, flattened or original
     if indexFlattenProcess > 0
         currentImg = double(imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, FileNames{1}{iFrame}]));
     else
         currentImg = double(movieData.channels_(iChannel).loadImage(iFrame));
     end
    
     level1 = thresholdOtsu(currentImg);
     currentImg = double(currentImg>level1).*currentImg;
     
    if iFrame==1
        % Get the first segmented results for the base of comparison
        current_img_firstframe = currentImg;
        int_sum_inside_firstframe = sum(sum(current_img_firstframe));
    end
    
    current_img_outside = double(currentImg).*(double(1-MaskFirstFrame));
    
    h12 = figure(12);
    hold off;
    
    imagesc(currentImg); colormap(gray);
    axis image; axis off;
    title(['Percentage of out growth based on intensity with zeroing the background: ', ...
        num2str(100*sum(sum(current_img_outside))/int_sum_inside_firstframe), '%'],'FontSize',20);
    
    hold on; plot(RoiYX(:,2),RoiYX(:,1),'m');
    
    saveas(h12,[HeatEnhOutputDir,'_bound/Int_display_',num2str(iFrame),'.tif']);
    
    int_outside_current_intensity(iChannel, iFrame) = sum(sum(current_img_outside));
    ratio_int_outside_firstframeinside(iChannel, iFrame) = sum(sum(current_img_outside))/int_sum_inside_firstframe;
    
end

display('The outgrowth percentage results based on intensity(with zeroing the background):');
display(ratio_int_outside_firstframeinside'*100);

% Save the outgrowth results
save([FilamentSegmentationChannelOutputDir,'/seg_outside.mat'], ...
    'seg_outside_current','seg_sum_inside_firstframe','ratio_outside_firstframeinside',...
    'ratio_int_outside_firstframeinside');
