function network_analysis(VIF_current_model,VIF_orientation, VIF_current_seg,outdir,iChannel,iFrame,ROI)
% function for calculation the property of network
% Liya Ding 01.2014.

% don't save every figure generated, unless debugging
save_everything_flag = 1;


img_size = size(VIF_current_seg);

VIF_current_seg = filament_model_to_seg_bwim(VIF_current_model,img_size,[]);

[VIF_digital_model,VIF_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
        = filament_model_to_digital_with_orientation(VIF_current_model);

    ROI = ones(img_size);
    ROI(:,1:300)=0;
    
    if(mean2(double(ROI))==1)
        VIF_ROI_model= VIF_digital_model;
        VIF_ROI_orientation_model =VIF_orientation_model;
    else
        count=1;
        VIF_ROI_model = cell(1,1);
        VIF_ROI_orientation_model =  cell(1,1);
        
        for iF = 1 : length(VIF_digital_model)
            x = VIF_digital_model{iF}(:,1);
            y = VIF_digital_model{iF}(:,2);
            ind_xy = sub2ind(img_size,y,x);
            roi_flag_array = ROI(ind_xy);
            [inROI_label, inROI_N] = bwlabel(roi_flag_array);
            for iR = 1 : inROI_N
               if (sum(double(inROI_label==iR))>2)
                   VIF_ROI_model{count} = [x(find(inROI_label==iR)) y(find(inROI_label==iR))];
                   VIF_ROI_orientation_model{count} = VIF_orientation_model{iF}(find(inROI_label==iR));                   
                   count = count +1 ;                   
               end
            end
        end
    end
    
%  [Vif_ROIed_model,Vif_ROIed_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
%             = filament_model_to_digital_with_orientation(VIF_ROI_model);
%        
    
orientation_pixel_pool = [];

pixel_number_per_filament_pool = [];
length_per_filament_pool = [];
straightness_per_filament_pool = [];



for iF = 1 : length(VIF_ROI_model)
    x = VIF_ROI_model{iF}(:,1);
    y = VIF_ROI_model{iF}(:,2);
    
    filament_detailed_length = sum(sqrt((x(1:end-1)-x(2:end)).^2 + (y(1:end-1)-y(2:end)).^2));
    filament_start_end_distance = sqrt((x(1)-x(end)).^2 + (y(1)-y(end)).^2);
    
    orientation_pixel_pool = [orientation_pixel_pool; VIF_ROI_orientation_model{iF}];
    pixel_number_per_filament_pool(iF) = length(x);
    length_per_filament_pool(iF) = filament_detailed_length;
    straightness_per_filament_pool(iF) = filament_start_end_distance/filament_detailed_length;
    
end

h1 =  figure(1); 

[h,bin] = hist(pixel_number_per_filament_pool,0:20:1000);
h = h/length(pixel_number_per_filament_pool);
bar(bin,h);
axis([-10 310 0 0.3]);


title('Pixels Distribution');
saveas(h1, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
saveas(h1, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
saveas(h1, [outdir,filesep,'network_pixels_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);

h2 =  figure(2); 


length_per_filament_pool = length_per_filament_pool(length_per_filament_pool>60);
[h,bin] = hist(length_per_filament_pool,0:20:1000);
h = h/length(length_per_filament_pool);
bar(bin,h);
axis([50 310 0 0.3]);


title('Length Distribution');
saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
saveas(h2, [outdir,filesep,'network_length_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);


orientation_pixel_pool_display = orientation_pixel_pool;

orientation_pixel_pool_display = orientation_pixel_pool_display + pi/2;
orientation_pixel_pool_display(orientation_pixel_pool_display>pi) = orientation_pixel_pool_display(orientation_pixel_pool_display>pi)-pi;
orientation_pixel_pool_display(orientation_pixel_pool_display>pi) = orientation_pixel_pool_display(orientation_pixel_pool_display>pi)-pi;
orientation_pixel_pool_display(orientation_pixel_pool_display>pi) = orientation_pixel_pool_display(orientation_pixel_pool_display>pi)-pi;
orientation_pixel_pool_display(orientation_pixel_pool_display<0) = orientation_pixel_pool_display(orientation_pixel_pool_display<0)+pi;
orientation_pixel_pool_display(orientation_pixel_pool_display<0) = orientation_pixel_pool_display(orientation_pixel_pool_display<0)+pi;
orientation_pixel_pool_display(orientation_pixel_pool_display<0) = orientation_pixel_pool_display(orientation_pixel_pool_display<0)+pi;

h3 =  figure(3); 
rose(orientation_pixel_pool_display);
title('Orientation of Filaments');
saveas(h3, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
saveas(h3, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
saveas(h3, [outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);

h6 =  figure(6); 

[h,bin] = hist(orientation_pixel_pool_display,0:pi/18:pi);
h = h/length(orientation_pixel_pool_display);
bar(bin,h);
axis([0-pi/36 pi+pi/36 0 0.3]);

title('Orientation of Filaments');
saveas(h6, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
saveas(h6, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
saveas(h6, [outdir,filesep,'network_orientationhist_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);



h4 =  figure(4); 
[h,bin] = hist(straightness_per_filament_pool,0:0.02:1);
h = h/length(straightness_per_filament_pool);
bar(bin,h);
axis([0.69 0.97 0 0.2]);

title('Straightness');
saveas(h4, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
saveas(h4, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
saveas(h4, [outdir,filesep,'network_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);

h5 =  figure(5); 
h = boxplot(straightness_per_filament_pool);
set(h(7,:),'Visible','off');
axis([0 2 0.6 1.01]);

title('Straightness');
saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.fig']);
saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.jpg']);
saveas(h5, [outdir,filesep,'network_box_straight_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'.tif']);


save([outdir,filesep,'network_orientationrose_ch_',num2str(iChannel),'_frame_',num2str(iFrame),'network_analysis.mat'], ...
    'straightness_per_filament_pool','orientation_pixel_pool_display',...
    'length_per_filament_pool');
