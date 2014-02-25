function load_MD_network_for_analysis(MD,ROI)

movie_Dir = MD.outputDirectory_;

package_process_ind_script;

for iChannel = 1 :  length(MD.channels_)
    
    outdir = [MD.processes_{indexFilamentSegmentationProcess}.outFilePaths_{iChannel},filesep,'analysis_results'];
    
    if(~exist(outdir,'dir'))
        mkdir(outdir);
    end
   
    for iFrame = 1 : nFrame
        display(['iChannel: ', num2str(iChannel),', iFrame:', num2str(iFrame)]);
        
        VIF_tip_orientation = MD.processes_{5}.loadChannelOutput(1,iFrame+0,'output','tip_orientation');
         VIF_RGB_seg_orient_heat_map = MD.processes_{5}.loadChannelOutput(1,iFrame+0,'output','RGB_seg_orient_heat_map');
        
         [vim_tip_y,vim_tip_x] = find(~isnan(VIF_tip_orientation));
         
         MT_tip_orientation = MD.processes_{5}.loadChannelOutput(2,iFrame+0,'output','tip_orientation');
         MT_RGB_seg_orient_heat_map = MD.processes_{5}.loadChannelOutput(2,iFrame+0,'output','RGB_seg_orient_heat_map');
        
          [mt_tip_y,mt_tip_x] = find(~isnan(MT_tip_orientation));
        
        CellSegmentation = MD.processes_{2}.loadChannelOutput(2,iFrame);
        
        CellBoundary = bwboundaries(CellSegmentation);
        CellBoundary = CellBoundary{1}; %% x in second col
        
           
        h1 = figure(1);
        imagesc(VIF_RGB_seg_orient_heat_map);axis off; axis image;
        hold on;
        plot(CellBoundary(:,2),CellBoundary(:,1),'r.');
        plot(vim_tip_x,vim_tip_y,'*');
        
        h2 = figure(2);
        imagesc(MT_RGB_seg_orient_heat_map);axis off; axis image;
        hold on;
        plot(CellBoundary(:,2),CellBoundary(:,1),'r.');
          plot(mt_tip_x,mt_tip_y,'r*');
      
        VIF_orientation = MD.processes_{5}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
        VIF_current_model = MD.processes_{5}.loadChannelOutput(iChannel,iFrame+0,'output','current_model');
        
        %     VIF_current_model
        [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
            = filament_model_to_digital_with_orientation(VIF_current_model);
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
%         ROI=[];
        if(isempty(ROI))
            ROI = ones(size(VIF_current_seg));
        end
        
%         network_analysis;       

        
        network_analysis(VIF_current_model,VIF_orientation, ...
                        VIF_current_seg,outdir,iChannel,iFrame,ROI);
        
        close all;
    end
end