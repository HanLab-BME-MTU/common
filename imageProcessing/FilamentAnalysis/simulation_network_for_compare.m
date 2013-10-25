movie_Dir = '/home/ld94/files/LCCB/interfil/Liya/MT_simulation_with_breaks/Processing/MT_release/4';

load([movie_Dir, filesep, 'movieData.mat']);

nFrame = MD.nFrames_;

SM_model = cell(1,nFrame);
MT_model = cell(1,nFrame);

flatten_dir{1} = MD.processes_{3}.outFilePaths_{1};
outdir = [MD.processes_{5}.outFilePaths_{1},filesep,'similarity_results'];
    mkdir(outdir);
    
    dist_pool_for_crossing=[];
    ang_pool_for_crossing=[];
    
for iFrame = 1 : 1
    iFrame
    
    % get the simulation ground truth with break lists
    load([movie_Dir, filesep,'MT_PSF_1.8_Vig_350_Noise_3_Breaks_2_BreakLength_2_3.54_Frame_100.mat']);
    afterrelease_MT_body_x_cell
    afterrelease_MT_body_y_cell
    break_list
    
    SM_current_model = cell(length(afterrelease_MT_body_x_cell),1);
    
    for iMT = 1 : length(afterrelease_MT_body_x_cell)
        SM_current_model{iMT,1}(:,1) = (afterrelease_MT_body_x_cell{iMT})';
        SM_current_model{iMT,1}(:,2) = (afterrelease_MT_body_y_cell{iMT})';
    end
%     SM_orientation = MD.processes_{5}.loadChannelOutput(1,iFrame+0,'output','current_seg_orientation');
%     SM_current_model = MD.processes_{5}.loadChannelOutput(1,iFrame+0,'output','current_model');
%     
%     [SM_digital_model,SM_orientation_model,SM_XX,SM_YY,SM_OO] ...
%     = filament_model_to_digital_with_orientation(SM_current_model);
% 
%      SM_current_seg = (isnan(SM_orientation)==0);
%      SM_img =  imread([flatten_dir{1},filesep,'flatten_',num2str(iFrame+0,'%03d'),'.tif']);
% 
%     % get the segmentation results
%     
%     MT_orientation = MD.processes_{5}.loadChannelOutput(1,iFrame,'output','current_seg_orientation');
% %     MT_orientation = MT_orientation(140:200,150:250);
    MT_current_model = MD.processes_{5}.loadChannelOutput(1,iFrame,'output','current_model');
%     
%     [MT_digital_model,MT_orientation_model,MT_XX,MT_YY,MT_OO] ...
%         = filament_model_to_digital_with_orientation(MT_current_model);
% 
%     
%     MT_current_seg = (isnan(MT_orientation)==0);
%     
    MT_img =  imread([flatten_dir{1},filesep,'flatten_f.tif']);
% %     MT_img = MT_img(140:200,150:250);
% %    SM_img = SM_img(140:200,150:250);
%    
%     
%     
%     two_channel_img = zeros(size(SM_img,1),size(SM_img,2),3);
%     two_channel_img(:,:,1)=SM_img;
%     two_channel_img(:,:,2)=MT_img;
%     h1=figure(1);imagesc(two_channel_img/255);axis equal;axis off;
%      saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.tif']);
%   saveas(h1,[outdir,filesep,'VIFMT_img_frame_',num2str(iFrame),'.fig']);
%    
    
    
%     two_channel_seg= zeros(size(SM_img,1),size(SM_img,2),3);
%     two_channel_seg(:,:,1)=SM_current_seg;
%      two_channel_seg(:,:,2)=MT_current_seg;
%     
%     h2=figure(2);imagesc(two_channel_seg);axis equal;axis off;
%    saveas(h2,[outdir,filesep,'VIFMT_seg_frame_',num2str(iFrame),'.tif']);
%    saveas(h2,[outdir,filesep,'VIFMT_seg_frame_',num2str(iFrame),'.fig']);
   
    network_similarity_with_breaks;
    
      close all;
end
