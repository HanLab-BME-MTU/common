function BA_output_cell = branch_analysis_movieList(ML,figure_flag)
% function to do branch analysis for a whole movielist
% Liya Ding, March, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function

% the number of movies
movieNumber =  length(ML.movieDataFile_);
BA_output_cell= cell(1,1);

for iM  = 1 :movieNumber
    
    clearvars -except 'movieNumber' 'BA_output_cell' 'iM' 'ML' 'figure_flag'
    
    % load this movie
    load(ML.movieDataFile_{iM});
    % the number of channels
    
    BA_output_cell = branch_analysis_movieData(MD,figure_flag);  
    
end

 ML_ROOT_DIR = ML.outputDirectory_;
 save([ML_ROOT_DIR,'\movieList_BA_output.mat'],'BA_output_cell');

 