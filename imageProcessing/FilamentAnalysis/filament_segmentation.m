function movieData = filament_segmentation(movieData, varargin)
% Created 07 2012 by Liya Ding, Matlab R2011b

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
    msg('Please set parameters for Image Flatten.')
    return;
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
Patch_Size = funParams.Patch_size;
Combine_Way = funParams.Combine_Way;

FilamentSegmentationOutputDir = funParams.OutputDirectory;

if (~exist(FilamentSegmentationOutputDir,'dir'))
    mkdir(FilamentSegmentationOutputDir);
end

nFrame = movieData.nFrames_;

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
    
    for iFrame = 1 : nFrame
        disp(['Frame: ',num2str(iFrame)]);
        
        % Read in the intensity image.
        if indexFlattenProcess > 0
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, FileNames{1}{iFrame}]);
         else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        
        load([SteerableChannelOutputDir, filesep, 'steerable_',num2str(iFrame),'.mat']);
        
        MaskCell = movieData.processes_{indexCellSegProcess}.loadChannelOutput(iChannel,iFrame);
       
        switch Combine_Way
            case 'int_st_both'
                
                 level0 = thresholdOtsu(MAX_st_res);
                 thresh_Segment = MAX_st_res > level0;
                                
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,Patch_Size,Pace_Size,0);
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,Patch_Size,Pace_Size,0);
                
                
            case 'st_only'                
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,Patch_Size,Pace_Size,0);
                
%                 [level1, SteerabelRes_Segment ] = thresholdOtsu_local_with_mask(MAX_st_res,MaskBig, Patch_Size,Pace_Size,0);
     
                current_seg = SteerabelRes_Segment; 
                Intensity_Segment = current_seg*0;
                
            case 'int_only'
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,Patch_Size,Pace_Size,0);
                % The segmentation is set as from intensity segmentation
                current_seg = Intensity_Segment; 
                SteerabelRes_Segment = current_seg*0;
            otherwise
                warning('Use the default of union');
                [level1, SteerabelRes_Segment ] = thresholdOtsu_local(MAX_st_res,Patch_Size,Pace_Size,0);
                [level2, Intensity_Segment ] = thresholdOtsu_local(currentImg,Patch_Size,Pace_Size,0);
                % The segmentation is set as the union of two segmentation.
                current_seg = or(Intensity_Segment,SteerabelRes_Segment);
        end
        

        current_seg = current_seg.*MaskCell;        
        
        current_seg_big = imdilate(current_seg, ones(5,5));
        current_seg_big = imfill(current_seg_big,'holes');
        
        %Label all objects in the mask
        labelMask = bwlabel(current_seg_big);
        
        %Get their area
        obAreas = regionprops(labelMask,'Area');      
        
        %First, check that there are objects to remove
        if length(obAreas) > 1
            obAreas = [obAreas.Area];
            %Sort by area
            [dummy,iSort] = sort(obAreas,'descend'); 
            %Keep only the largest requested number
            current_seg_big = labelMask == iSort(1);
        end
        
        current_seg_big = imerode(current_seg_big,ones(5,5));        
        current_seg = current_seg_big.*current_seg;
        
        
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

%         imwrite(current_seg,[FilamentSegmentationChannelOutputDir,'/segment_',num2str(iFrame),'.tif']);
        
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
        
        save([FilamentSegmentationChannelOutputDir,'/steerable_vote_',num2str(iFrame),'.mat'],...
            'currentImg','orienation_map_filtered','OrientationVoted','orienation_map', ...
            'MAX_st_res', 'current_seg','Intensity_Segment','SteerabelRes_Segment');
        
    end
end     
 