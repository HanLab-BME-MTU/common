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
        
        VIF_orientation = MD.processes_{5}.loadChannelOutput(iChannel,iFrame+0,'output','current_seg_orientation');
        VIF_current_model = MD.processes_{5}.loadChannelOutput(iChannel,iFrame+0,'output','current_model');
        
        %     VIF_current_model
        [Vif_digital_model,Vif_orientation_model,VIF_XX,VIF_YY,VIF_OO] ...
            = filament_model_to_digital_with_orientation(VIF_current_model);
        
        VIF_current_seg = (isnan(VIF_orientation)==0);
        ROI=[];
        if(isempty(ROI))
            ROI = ones(size(VIF_current_seg));
        end
        
%         network_analysis;       

        
        network_analysis(VIF_current_model,VIF_orientation, ...
                        VIF_current_seg,outdir,iChannel,iFrame,ROI);
        
        close all;
    end
end