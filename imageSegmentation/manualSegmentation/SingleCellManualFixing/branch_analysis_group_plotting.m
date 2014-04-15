function BA_Group = branch_analysis_group_plotting(ML_name_cell, T_branchsize, T_branchduration, Group_ROOT_DIR)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML_array:     The array of movieList objects in this group to be analyzed

nList = length(ML_name_cell);

BA_Group = cell(1, nList);

if(T_branchsize==0)
    T_branchsize=Inf;
end

if(T_branchduration==0)
    T_branchduration=Inf;
end

close all;

for iML = 1 : nList
    
    ML_name = ML_name_cell{iML};
    
    load(ML_name);
    
    ML_ROOT_DIR = ML.outputDirectory_;
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    % if the batch results exist, load it
%     if(exist([ML_ROOT_DIR,'\movieList_BA_output.mat'], 'file'))
    if(0)
        load([ML_ROOT_DIR,'\movieList_BA_output.mat'],'BA_output_ML_cell');
    else
        % otherwise load one by one
        BA_output_ML_cell= cell(1,1,1);
        
        for iM  = 1 :movieNumber
            % load this movie
            load(ML.movieDataFile_{iM});
            
            % Now MD in workspace
            ROOT_DIR = MD.outputDirectory_;
            
            % the number of channels
            nChannel = length(MD.channels_);
            
            for iChannel = 1 : nChannel
                for iCell = 1 : 20
                    
                    display(['iML: ',num2str(iML), ', iM:', num2str(iM), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
                    
                    
                    % the folder name if there is marking
                    outputPath = [ROOT_DIR,'\BranchAnalysisChannel',num2str(iChannel),'Cell',num2str(iCell)];
                    
                    % check if this folder exist
                    if(exist(outputPath,'dir'))
                        % if it exist, try to do the branch analysis
                        try
                            
                            if(exist( [outputPath,'\branch_analysis_results.mat'],'file'))
                                load([outputPath,'\branch_analysis_results.mat'],'BA_output');
                            else
                                %                                 BA_output = branch_analysis_marked_cell(MD, iChannel, iCell);
                                continue;
                            end
                            
                            % for the output, ignore the channel number since
                            % there is only one channel marked.
                            BA_output_ML_cell{1, iM}{iChannel, iCell} = BA_output;
                        end
                    end
                end
            end
        end
    end
    
    BA_Group{1, iML} = BA_output_ML_cell;
    
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
Group_Pool_whole_cell_size_mean= [];
Group_Pool_branch_number_mean_pat=[];
Group_Pool_whole_cell_vim_totalamount_mean=[];
Group_Pool_whole_cell_vif_mean_intensity_pat=[];
Group_Pool_branch_size_mean=[];

iAllCell = 0;

if(isempty(Group_ROOT_DIR))
    Group_ROOT_DIR = ML_ROOT_DIR;
end

fid = fopen([Group_ROOT_DIR,'\branch_analysis_group_list.txt'],'wt');

fclose(fid);

for iML = 1 : nList
    
    ML_name = ML_name_cell{iML};
    load(ML_name);
    ML_ROOT_DIR = ML.outputDirectory_;
    
    BA_output_ML_cell = BA_Group{1, iML};
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    for iM  = 1 : movieNumber
        
        % load this movie
        load(ML.movieDataFile_{iM});
        MD_dir = MD.outputDirectory_;
        
        % the number of channels
        nChannel = length(MD.channels_);
        
        for iChannel = 1 : nChannel
            try
                for iCell = 1 :  size(BA_output_ML_cell{1, iM},2)
                    try
                        BA_output = BA_output_ML_cell{1, iM}{iChannel, iCell};
                        if(~isempty(BA_output))
                            
                            
                            Group_Pool_branch_size_mean = [Group_Pool_branch_size_mean BA_output.branch_mean_size];
                            Group_Pool_whole_cell_size_mean= [Group_Pool_whole_cell_size_mean BA_output.whole_cell_size_mean];
                            Group_Pool_branch_number_mean_pat = [Group_Pool_branch_number_mean_pat repmat(BA_output.branch_number_mean, [1, length(BA_output.branch_vif_mean_intensity)])];
                            Group_Pool_whole_cell_vif_mean_intensity_pat = [Group_Pool_whole_cell_vif_mean_intensity_pat repmat(BA_output.whole_cell_vif_mean_intensity, [1, length(BA_output.branch_vif_mean_intensity)])];
                            Group_Pool_whole_cell_vim_totalamount_mean = [Group_Pool_whole_cell_vim_totalamount_mean BA_output.whole_cell_vim_totalamount_mean];
                            
                            
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
                            
                            iAllCell = iAllCell + 1;
                            
                            if(MD_dir(end)=='/' || MD_dir(end)=='\')
                                MD_dir = MD_dir(1:end-1);
                            end
                                                       
                            
                            ind = find(MD_dir(1:end-1)=='/' | MD_dir(1:end-1)=='\' );
                            
                            
                            movie_dir_list{iAllCell} = [MD_dir(ind(end)+1:end), '     Channel: ', num2str(iChannel),',  Cell: ', num2str(iCell)];
                            
                            str_line = ['FinalCell Index: ', num2str(iAllCell),':  iML: ',num2str(iML), ', iM:', num2str(iM),', ', MD_dir(ind(end)+1:end), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)];
                            display(['FinalCell Index: ', num2str(iAllCell),':  iML: ',num2str(iML), ', iM:', num2str(iM),', ', MD_dir(ind(end)+1:end), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
                            fprintf(fid, [str_line,'\n  \r\n']);


                        end
                    end
                end
            end
        end
    end
end
Group_Pool_tracked_branch_d_frames = Group_Pool_branch_number_tracked./Group_Pool_Cell_Marked_Frame_Number;

colorarray = rand(length(Group_Pool_branch_vif_mean_intensity),3);

%% % plot the this group result
h3 = figure(3);
plot(Group_Pool_branch_duration_array, Group_Pool_branch_vif_mean_intensity,'.');
xlabel('Branch Duration');
ylabel('Branch Vim mean Int');
title('Branch: Duration vs Vim');
saveas(h3,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim.fig']);
saveas(h3,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim.jpg']);

h4 = figure(4);
[h_p,bin] = hist(Group_Pool_protrusion_vif_mean_intensity,100:20:300);
[h_r,bin] = hist(Group_Pool_retraction_vif_mean_intensity,100:20:300);
bar([reshape(bin, 1, numel(bin)); reshape(bin, 1, numel(bin))]',[reshape(h_p, 1, numel(h_p));reshape(h_r, 1, numel(h_r))]');
legend('Protrusion Region Vim','Retraction Region Vim');
saveas(h4,[Group_ROOT_DIR,'\Protrusion_vs_Retraction.fig']);
saveas(h4,[Group_ROOT_DIR,'\Protrusion_vs_Retraction.jpg']);

% mean vim vs speed
h7 = figure(7); hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_Travel_Speed,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_Travel_Speed(i)+0.2, num2str(i), 'color',colorarray(i,:))
hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), '.', 'color',colorarray(i,:))
end
xlabel('Branch Mean Number');
ylabel('Cell Speed');
title(['Branch: Branchness vs Speed (Sample number: ', num2str(length(Group_Pool_whole_cell_vif_mean_intensity)),')']);
saveas(h7,[Group_ROOT_DIR,'\Branchness_vs_Speed.fig']);
saveas(h7,[Group_ROOT_DIR,'\Branchness_vs_Speed.jpg']);

h8 = figure(8);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
hold on;
plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:))
end
xlabel('Branch Mean Number');
ylabel('Cell Vim mean Int');
title(['Branch number vs Vim (Sample number: ', num2str(length(Group_Pool_whole_cell_vif_mean_intensity)),')']);
saveas(h8,[Group_ROOT_DIR,'\Branchness_vs_Vim.fig']);
saveas(h8,[Group_ROOT_DIR,'\Branchness_vs_Vim.jpg']);



h9 = figure(9);
plot(Group_Pool_branch_number_mean_pat, Group_Pool_branch_vif_mean_intensity,'.');
xlabel('Branch Mean Number');
ylabel('Branch Vim mean Int');
title('Branch number vs Vim');
saveas(h9,[Group_ROOT_DIR,'\Branchness_vs_BranchVim.fig']);
saveas(h9,[Group_ROOT_DIR,'\Branchness_vs_BranchVim.jpg']);





h10 = figure(10);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vim_totalamount_mean)
text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
hold on;
plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
end
xlabel('Branch Mean Number');
ylabel('Cell Vim Total Int');
title(['Branch number vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')']);
saveas(h10,[Group_ROOT_DIR,'\Branchness_vs_VimTotal.fig']);
saveas(h10,[Group_ROOT_DIR,'\Branchness_vs_VimTotal.jpg']);



h11 = figure(11);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_tracked_branch_d_frames)
text(Group_Pool_tracked_branch_d_frames(i)-0.02, Group_Pool_whole_cell_vim_totalamount_mean(i)+10, num2str(i), 'color',colorarray(i,:));
hold on;
plot(Group_Pool_tracked_branch_d_frames(i), Group_Pool_whole_cell_vim_totalamount_mean(i), '.', 'color',colorarray(i,:))
end
xlabel('Tracked Branchs per Frame');
ylabel('Cell Vim Total Int');
title(['Tracked Branches per Frame vs Vim Total(Sample number: ', num2str(length(Group_Pool_whole_cell_vim_totalamount_mean)),')']);
saveas(h11,[Group_ROOT_DIR,'\TrackedBranchness_vs_VimTotal.fig']);
saveas(h11,[Group_ROOT_DIR,'\TrackedBranchness_vs_VimTotal.jpg']);



h12 = figure(12);
plot(Group_Pool_whole_cell_vif_mean_intensity_pat, Group_Pool_branch_duration_array, '.');
ylabel('Branch Duration');
xlabel('Cell Vim mean Int');
title('Cell Vim vs Branch Duration');
saveas(h12,[Group_ROOT_DIR,'\BranchDuration_vs_CellVim.fig']);
saveas(h12,[Group_ROOT_DIR,'\BranchDuration_vs_CellVim.jpg']);

h13 = figure(13);
plot(Group_Pool_branch_duration_array(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)), ...
    Group_Pool_branch_vif_mean_intensity(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)),'.');
xlabel('Branch Duration');
ylabel('Branch Vim mean Int');
title(['Branch: Duration vs Vim, for branches larger than ',num2str(T_branchsize),' pixels, duration longer than ',num2str(T_branchduration),' frames']);
saveas(h13,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim_with_Thres.fig']);
saveas(h13,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim_with_Thres.jpg']);



