% function network_feature_ML_cell = load_ML_network_analysis(ML,figure_flag)
% function to do network analysis for a whole movielist
% Liya Ding, June, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function

% the number of movies

movieNumber =  length(ML.movieDataFile_);
network_feature_ML_cell= cell(1,1);

for iM  = 1 :movieNumber
    
    clearvars -except 'movieNumber' 'network_feature_ML_cell' 'iM' 'ML' 'figure_flag'
    
    close all;
    
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    display('======================================================================');
    display(['iM:', num2str(iM)]);
                     
    network_feature_ML_cell{1, iM} = load_MD_network_for_analysis(MD,[]);  
    
end

 ML_ROOT_DIR = ML.outputDirectory_;
 save([ML_ROOT_DIR,'\movieList_netwrok_analysis_output.mat'],'network_feature_ML_cell');

 
 %% plot the density out
 
 MT_density_pool = [];
 VIF_density_pool = [];
 
 
for iM  = 1 :movieNumber
    this_MD_f = network_feature_ML_cell{1, iM};
    nFrame = size(this_MD_f,2);
    
%     for iF = 1 : nFrame
%         MT_density = this_MD_f{1,2}.density_filament;
%         VIF_density = this_MD_f{1,1}.density_filament;
%         
%         MT_cell_seg = this_MD_f{1,2}.Cell_Mask;
%         VIF_cell_seg = this_MD_f{1,1}.Cell_Mask;
%         Cell_Mask = MT_cell_seg.*VIF_cell_seg;
%         
%         MT_density_pool = [MT_density_pool MT_density(Cell_Mask>0)];
%         VIF_density_pool = [VIF_density_pool VIF_density(Cell_Mask>0)];   
%     end    
%     
    
     for iF = 1 : nFrame
        MT_density = this_MD_f{1,iF}.density_filament;
        VIF_density = this_MD_f{2,iF}.density_filament;
        
        MT_cell_seg = this_MD_f{1,iF}.Cell_Mask;
        VIF_cell_seg = this_MD_f{2,iF}.Cell_Mask;
        
        
        Cell_Mask = MT_cell_seg.*VIF_cell_seg.*(MT_density>0.002).*(VIF_density>0.002);
        
        MT_density_pool = [MT_density_pool; MT_density(Cell_Mask>0)];
        VIF_density_pool = [VIF_density_pool; VIF_density(Cell_Mask>0)];   
    end   
end

figure;
plot(MT_density_pool(1:100:end),VIF_density_pool(1:100:end),'.');
title(['MT and VIF densities Cross-correlation:', num2str(corr(MT_density_pool,VIF_density_pool))]);

 