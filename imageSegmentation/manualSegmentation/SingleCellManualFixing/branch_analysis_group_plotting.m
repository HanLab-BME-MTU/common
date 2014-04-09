% function BA_output_cell = branch_analysis_group_plotting(ML_name_cell)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML_array:     The array of movieList objects in this group to be analyzed

nList = length(ML_name_cell);

BA_Group = cell(1, nList);

for iML = 1 : nList
    
    ML_name = ML_name_cell{iML};
    
    load(ML_name);
    
    ML_ROOT_DIR = ML.outputDirectory_;
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    % if the batch results exist, load it
    if(exist([ML_ROOT_DIR,'\movieList_BA_output.mat'], 'file'))
        load([ML_ROOT_DIR,'\movieList_BA_output.mat'],'BA_output_cell');
    else
        % otherwise load one by one   
        BA_output_cell= cell(1,1);
        
        for iM  = 1 :movieNumber
            % load this movie
            load(ML.movieDataFile_{iM});
            
            % Now MD in workspace
            ROOT_DIR = MD.outputDirectory_;
            
            % the number of channels
            nChannel = length(MD.channels_);
            
            for iChannel = 1 : nChannel
                for iCell = 1 : 10
                    % the folder name if there is marking
                    outputPath = [ROOT_DIR,'\BranchAnalysisChannel',num2str(iChannel),'Cell',num2str(iCell)];
                    
                    % check if this folder exist
                    if(exist(outputPath,'dir'))
                        % if it exist, try to do the branch analysis
                        try
                            
                            if(exist( [outputPath,'\branch_analysis_results.mat'],'file'))
                                load([outputPath,'\branch_analysis_results.mat'],'BA_output');
                            else
                                BA_output = branch_analysis_marked_cell(MD, iChannel, iCell);
                            end
                            
                            % for the output, ignore the channel number since
                            % there is only one channel marked.
                            BA_output_cell{iM, iCell} = BA_output;
                        end
                    end
                end
            end
        end
    end
    
    BA_Group{1, iML} = BA_output_cell;

end

% initialize the pools

Group_Pool_Travel_Length = [];
Group_Pool_Travel_Distance = [];
Group_Pool_Travel_Speed = [];
Group_Pool_Cell_Marked_Frame_Number = [];
Group_Pool_branch_number_tracked = [];
Group_Pool_branch_number_max = [];
Group_Pool_branch_number_mean = [];
Group_Pool_branch_duration_array = [];
Group_Pool_branch_duration_mean= [];
Group_Pool_branch_vif_mean_intensity= [];
Group_Pool_protrusion_vif_mean_intensity= [];
Group_Pool_retraction_vif_mean_intensity= [];
Group_Pool_whole_cell_vif_mean_intensity= [];


for iML = 1 : nList
     
    ML_name = ML_name_cell{iML};    
    load(ML_name);
      ML_ROOT_DIR = ML.outputDirectory_;
    
    BA_output_cell = BA_Group{1, iML};
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    for iM  = 1 : movieNumber
        for iCell = 1 :  size(BA_output_cell,2)
             BA_output = BA_output_cell{iM, iCell};
             if(~isempty(BA_output))
                 
                 
                 Group_Pool_Travel_Length = [Group_Pool_Travel_Length BA_output.cell_travel_length];
                 Group_Pool_Travel_Distance = [Group_Pool_Travel_Distance BA_output.cell_travel_distance];
                 Group_Pool_Travel_Speed = [Group_Pool_Travel_Speed BA_output.cell_travel_speed];
                 Group_Pool_Cell_Marked_Frame_Number = [Group_Pool_Cell_Marked_Frame_Number BA_output.cell_marked_frame_number];
                 Group_Pool_branch_number_tracked = [Group_Pool_branch_number_tracked BA_output.branch_number_tracked];
                 Group_Pool_branch_number_max = [Group_Pool_branch_number_max BA_output.branch_number_max];
                 Group_Pool_branch_number_mean = [Group_Pool_branch_number_mean BA_output.branch_number_mean];
                 Group_Pool_branch_duration_array = [Group_Pool_branch_duration_array  BA_output.branch_duration_array];
                 Group_Pool_branch_duration_mean= [Group_Pool_branch_duration_mean  BA_output.branch_duration_mean];
                 Group_Pool_branch_vif_mean_intensity= [Group_Pool_branch_vif_mean_intensity BA_output.branch_vif_mean_intensity];
                 Group_Pool_protrusion_vif_mean_intensity= [Group_Pool_protrusion_vif_mean_intensity BA_output.protrusion_vif_mean_intensity];
                 Group_Pool_retraction_vif_mean_intensity= [Group_Pool_retraction_vif_mean_intensity BA_output.retraction_vif_mean_intensity];
                 Group_Pool_whole_cell_vif_mean_intensity= [Group_Pool_whole_cell_vif_mean_intensity BA_output.whole_cell_vif_mean_intensity];

             end
        end
    end                          
end

%% % plot the this group result
figure(3);
plot(Group_Pool_branch_duration_array, Group_Pool_branch_vif_mean_intensity,'.');
xlabel('Branch Duration');
ylabel('Branch Vim mean Int');
title('Branch: Duration vs Vim');

figure(4);
[h_p,bin] = hist(Group_Pool_protrusion_vif_mean_intensity,100:50:300);
[h_r,bin] = hist(Group_Pool_retraction_vif_mean_intensity,100:50:300);
bar([bin; bin]',[h_p;h_r]');
legend('Protrusion Region Vim','Retraction Region Vim');

figure(5);
plot(Group_Pool_Travel_Speed, Group_Pool_branch_number_tracked,'.');
xlabel('Cell Speed');
ylabel('Branch numbers');
title('Motility with Branchiness');

