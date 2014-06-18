function BA_Group = branch_analysis_group_plotting(ML_name_cell, T_branchsize, T_branchduration, Group_ROOT_DIR)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML_array:     The array of movieList objects in this group to be analyzed

nList = length(ML_name_cell);

BA_Group = cell(1, nList);

if(nargin<2)
    T_branchsize=200;
end
if(nargin<3)
    T_branchduration=2;
end
if(nargin<4)
    Group_ROOT_DIR=[];
end

close all;

for iML = 1 : nList
    
    ML_name = ML_name_cell{iML};
    
    load(ML_name);
    
    ML_ROOT_DIR = ML.outputDirectory_;
    
    % the number of movies
    movieNumber =  length(ML.movieDataFile_);
    
    % if the batch results exist, load it
    if(exist([ML_ROOT_DIR,'\movieList_BA_results_gathered.mat'], 'file'))
        load([ML_ROOT_DIR,'\movieList_BA_results_gathered.mat'],'BA_output_ML_cell');
    else
        % otherwise load one by one
        BA_output_ML_cell= cell(1,1,1);
        
        for iM  = 1 :movieNumber
            try
                % load this movie
                load(ML.movieDataFile_{iM});
            catch
                disp('The MD file is missing');
                continue;
            end
            % Now MD in workspace
            ROOT_DIR = MD.outputDirectory_;
            
            % the number of channels
            nChannel = length(MD.channels_);
            
            for iChannel = 1 : nChannel
                for iCell = 1 : 20
                    display(['Checking: iM:', num2str(iM), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
                    
                    % the folder name if there is marking
                    outputPath = [ROOT_DIR,'\BranchAnalysisChannel',num2str(iChannel),'Cell',num2str(iCell)];
                    
                    % check if this folder exist
                    if(exist(outputPath,'dir'))
                        % if it exist, try to do the branch analysis
                        try
                            if(exist([outputPath,'\branch_analysis_results.mat'],'file'))
                                load([outputPath,'\branch_analysis_results.mat'],'BA_output');
                            else
                                continue;
                            end
                            
                            display(['Found: iM:', num2str(iM), ', iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
                            
                            % save
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
Group_Pool_thresholded_branch_number_mean=[];
Group_Pool_branch_cellmovement_std=[];

iAllCell = 0;

% if no input, the DIR of last movieList
if(isempty(Group_ROOT_DIR))
    Group_ROOT_DIR = ML_ROOT_DIR;
end

if(exist(Group_ROOT_DIR,'dir')==0)
    mkdir(Group_ROOT_DIR);
end

for iML = 1 : nList
    BA_output_ML_cell = BA_Group{1, iML};
    
    % the number of movies
    movieNumber =  20;
    
    for iM  = 1 : movieNumber
        nChannel=2;         
        for iChannel = 1 : nChannel
            for iCell = 1 :  10
                try
                    BA_output = BA_output_ML_cell{1, iM}{iChannel, iCell};
                catch
                    continue;
                end
                
                if(~isempty(BA_output))                            
                            Group_Pool_branch_size_mean = [Group_Pool_branch_size_mean BA_output.branch_mean_size];
                            Group_Pool_whole_cell_size_mean= [Group_Pool_whole_cell_size_mean BA_output.whole_cell_size_mean];
                            Group_Pool_branch_number_mean_pat = [Group_Pool_branch_number_mean_pat repmat(BA_output.branch_number_mean, [1, length(BA_output.branch_vif_mean_intensity)])];
                            Group_Pool_whole_cell_vif_mean_intensity_pat = [Group_Pool_whole_cell_vif_mean_intensity_pat repmat(BA_output.whole_cell_vif_mean_intensity, [1, length(BA_output.branch_vif_mean_intensity)])];
                            Group_Pool_whole_cell_vim_totalamount_mean = [Group_Pool_whole_cell_vim_totalamount_mean BA_output.whole_cell_vim_totalamount_mean];
                            
                            
                            Thresholded_branch_number_mean = sum(BA_output.branch_duration_array(find(BA_output.branch_mean_size>T_branchsize & BA_output.branch_duration_array>T_branchduration)))...
                                        ./ BA_output.cell_marked_frame_number;
                            Group_Pool_thresholded_branch_number_mean = [Group_Pool_thresholded_branch_number_mean Thresholded_branch_number_mean];
                                                     
                            
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
                            Group_Pool_branch_cellmovement_std= [Group_Pool_branch_cellmovement_std BA_output.branch_cellmovement_std];
                           
                            iAllCell = iAllCell + 1;
                            
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
xlabel('Branch Duration','Fontsize',14);
ylabel('Branch Vim mean Int','Fontsize',14);
title('Branch: Duration vs Vim','Fontsize',14);
saveas(h3,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim.fig']);
saveas(h3,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim.tif']);

h4 = figure(4);
[h_p,bin] = hist(Group_Pool_protrusion_vif_mean_intensity,100:20:500);
[h_r,bin] = hist(Group_Pool_retraction_vif_mean_intensity,100:20:500);
bar([reshape(bin, 1, numel(bin)); reshape(bin, 1, numel(bin))]',[reshape(h_p, 1, numel(h_p));reshape(h_r, 1, numel(h_r))]');
legend('Protrusion Region Vim','Retraction Region Vim');
saveas(h4,[Group_ROOT_DIR,'\Protrusion_vs_Retraction.fig']);
saveas(h4,[Group_ROOT_DIR,'\Protrusion_vs_Retraction.tif']);

% mean vim vs speed
h7 = figure(7); hold off;
%  plot(Group_Pool_branch_number_mean, Group_Pool_Travel_Speed,'o');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_Travel_Speed(i)+0.2, num2str(i), 'color',colorarray(i,:))
% hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'o', 'color',colorarray(i,:))
hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'bo', 'MarkerSize',6);
hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'bo', 'MarkerSize',7);
hold on;plot(Group_Pool_branch_number_mean(i), Group_Pool_Travel_Speed(i), 'bo', 'MarkerSize',8);
end
xlabel('Average Branch Number','FontSize',14);
ylabel('Cell Speed','FontSize',14);
title(['Branchness vs Speed, ',...
        ' correlation: ',...
    num2str(corr(Group_Pool_branch_number_mean', ...
    Group_Pool_Travel_Speed'),'%1.2f')],'FontSize',14);
% axis([0 10 0 24]);
saveas(h7,[Group_ROOT_DIR,'\Branchness_vs_Speed.fig']);
saveas(h7,[Group_ROOT_DIR,'\Branchness_vs_Speed.tif']);
saveas(h7,[Group_ROOT_DIR,'\Branchness_vs_Speed.tif']);

h8 = figure(8);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
hold on;
plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16)
end
xlabel('Branch Mean Number','Fontsize',14);
ylabel('Cell Vim mean Int','Fontsize',14);
title(['Branch number vs Vim, correlation: ',...
    num2str(corr(Group_Pool_branch_number_mean', ...
    Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',14);
saveas(h8,[Group_ROOT_DIR,'\Branchness_vs_Vim.fig']);
saveas(h8,[Group_ROOT_DIR,'\Branchness_vs_Vim.tif']);
h8_axis = axis;

h18 = figure(18);hold off;
% plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
text(Group_Pool_thresholded_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
hold on;
plot(Group_Pool_thresholded_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
end
xlabel('Branch Mean Number','Fontsize',14);
ylabel('Cell Vim mean Int','Fontsize',14);
title(['Branch number vs Vim, for branches larger than ',num2str(T_branchsize),' pixels, duration longer than ',num2str(T_branchduration),' frames']);
axis(h8_axis);
saveas(h18,[Group_ROOT_DIR,'\Branchness_vs_Vim_Thresholded.fig']);
saveas(h18,[Group_ROOT_DIR,'\Branchness_vs_Vim_Thresholded.tif']);
corr_b_v_v = corrcoef(Group_Pool_whole_cell_vif_mean_intensity',Group_Pool_thresholded_branch_number_mean');
display(['Corr for Branchness vs Vim: ',num2str(corr_b_v_v(1,2))]);

%%
h28 = figure(28);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
hold on;
% plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','MarkerSize',6);
plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','MarkerSize',7);
plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_whole_cell_vif_mean_intensity(i), 'bo','MarkerSize',8);
end
% axis([0.55, 1.0, 350 950]);
xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',14);
ylabel('Cell Vim mean Int','Fontsize',14);
title(['Branch Orientation Scatterness vs Vim ,',...
    ' correlation: ',...
    num2str(corr(Group_Pool_branch_cellmovement_std', ...
    Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',14);
saveas(h28,[Group_ROOT_DIR,'\BranchOrient_vs_Vim.fig']);
saveas(h28,[Group_ROOT_DIR,'\BranchOrient_vs_Vim.tif']);
saveas(h28,[Group_ROOT_DIR,'\BranchOrient_vs_Vim.tif']);
%%

h38 = figure(38);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_branch_number_mean)
% text(Group_Pool_branch_cellmovement_std(i)-0.02, Group_Pool_branch_number_mean(i)+0.2, num2str(i), 'color',colorarray(i,:));
hold on;
plot(Group_Pool_branch_cellmovement_std(i), Group_Pool_branch_number_mean(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
end
xlabel('Branch Orientation along Cell Movement, Standard Deviation','Fontsize',14);
ylabel('Branch Mean Number ','Fontsize',14);
title(['Branch Orientation Scatterness vs Branch Mean Number,',...
    ' correlation: ',...
    num2str(corr(Group_Pool_branch_cellmovement_std', ...
    Group_Pool_branch_number_mean'),'%1.2f') ],'Fontsize',14);

saveas(h38,[Group_ROOT_DIR,'\BranchOrient_vs_Branchness.fig']);
saveas(h38,[Group_ROOT_DIR,'\BranchOrient_vs_Branchness.tif']);

h48 = figure(48);hold off;
% plot(Group_Pool_branch_cellmovement_std, Group_Pool_whole_cell_vif_mean_intensity,'.');
for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% text(Group_Pool_Travel_Speed(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+0.2, num2str(i), 'color',colorarray(i,:));
hold on;
plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',16);
end
xlabel('Cell Speed','Fontsize',14);
ylabel('Vim Level ','Fontsize',14);
title(['Cell Speed vs Vim Level, ',...
    ' correlation: ',...
    num2str(corr(Group_Pool_Travel_Speed', ...
    Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',14);

saveas(h48,[Group_ROOT_DIR,'\Speed_vs_Vim.fig']);
saveas(h48,[Group_ROOT_DIR,'\Speed_vs_Vim.tif']);

% 
h9 = figure(9);
plot(Group_Pool_branch_number_mean_pat, Group_Pool_branch_vif_mean_intensity,'.');
xlabel('Branch Mean Number');
ylabel('Branch Vim mean Int');
title('Branch number vs Vim');
saveas(h9,[Group_ROOT_DIR,'\Branchness_vs_BranchVim.fig']);
saveas(h9,[Group_ROOT_DIR,'\Branchness_vs_BranchVim.tif']);





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
saveas(h10,[Group_ROOT_DIR,'\Branchness_vs_VimTotal.tif']);


% h20 = figure(20);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% hold on;
% plot(Group_Pool_branch_number_mean(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% end
% xlabel('Branch Mean Number','Fontsize',14);
% ylabel('Cell Vim Mean Int','Fontsize',14);
% title(['Branch number vs Vim Mean,  correlation: ',...
%     num2str(corr(Group_Pool_branch_number_mean', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',14);
% saveas(h20,[Group_ROOT_DIR,'\Branchness_vs_VimMean.fig']);
% saveas(h20,[Group_ROOT_DIR,'\Branchness_vs_VimMean.tif']);
% 
% h30 = figure(30);hold off;
% % plot(Group_Pool_branch_number_mean, Group_Pool_whole_cell_vif_mean_intensity,'.');
% for i = 1 : length(Group_Pool_whole_cell_vif_mean_intensity)
% % text(Group_Pool_branch_number_mean(i)-0.02, Group_Pool_whole_cell_vif_mean_intensity(i)+10, num2str(i), 'color',colorarray(i,:));
% hold on;
% plot(Group_Pool_Travel_Speed(i), Group_Pool_whole_cell_vif_mean_intensity(i), '.', 'color',colorarray(i,:),'MarkerSize',15)
% end
% xlabel('Cell Speed','Fontsize',14);
% ylabel('Cell Vim Mean Int','Fontsize',14);
% title(['Cell Speed vs Vim Mean,  correlation: ',...
%     num2str(corr(Group_Pool_Travel_Speed', ...
%     Group_Pool_whole_cell_vif_mean_intensity'),'%1.2f') ],'Fontsize',14);
% saveas(h30,[Group_ROOT_DIR,'\Speed_vs_VimMean.fig']);
% saveas(h30,[Group_ROOT_DIR,'\Speed_vs_VimMean.tif']);


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
saveas(h11,[Group_ROOT_DIR,'\TrackedBranchness_vs_VimTotal.tif']);


h12 = figure(12);
plot(Group_Pool_whole_cell_vif_mean_intensity_pat, Group_Pool_branch_duration_array, '.');
ylabel('Branch Duration');
xlabel('Cell Vim mean Int');
title('Cell Vim vs Branch Duration');
saveas(h12,[Group_ROOT_DIR,'\BranchDuration_vs_CellVim.fig']);
saveas(h12,[Group_ROOT_DIR,'\BranchDuration_vs_CellVim.tif']);

h13 = figure(13);
plot(Group_Pool_branch_duration_array(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)), ...
    Group_Pool_branch_vif_mean_intensity(find(Group_Pool_branch_size_mean>T_branchsize & Group_Pool_branch_duration_array>T_branchduration)),'.');
xlabel('Branch Duration');
ylabel('Branch Vim mean Int');
title(['Branch: Duration vs Vim, for branches larger than ',num2str(T_branchsize),' pixels, duration longer than ',num2str(T_branchduration),' frames']);
saveas(h13,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim_with_Thres.fig']);
saveas(h13,[Group_ROOT_DIR,'\EachBranch_Duration_vs_Vim_with_Thres.tif']);


save([Group_ROOT_DIR,'\results.mat']);
