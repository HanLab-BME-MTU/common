function updateMovieData(movieData)

% updateMovieData(movieData)
%
% This will save the input movieData to it's analysis directory, and backup
% the previous movieData as movieData.old
%
% The input movieData can be a single movieData structure or a cell-array
% of movieDatas (a.k.a. a movieArray)
%
% Hunter Elliott


%Convert to cell if a structure was input
if ~iscell(movieData)
    movieArray = {movieData};
else
    movieArray = movieData;
end

%Get number of movies
nMovies = length(movieArray);

for j = 1:nMovies
    if isfield(movieArray{j},'analysisDirectory') && ...
            exist(movieArray{j}.analysisDirectory,'dir')

        if exist([movieArray{j}.analysisDirectory filesep 'movieData.mat'],'file') %Keep the old one
            
             copyfile([movieArray{j}.analysisDirectory filesep 'movieData.mat'], ...
                      [movieArray{j}.analysisDirectory filesep 'movieData.old']);
            
        end
        
        movieData = movieArray{j}; %#ok<NASGU>
        save([movieArray{j}.analysisDirectory filesep 'movieData.mat'],'movieData');
        if exist('oldMovieData','var')            
            movieData = oldMovieData;
            save([movieArray{j}.analysisDirectory filesep 'movieData.old'],'movieData','-mat');
        end

    else
        error('Cannot locate movie analysis directory!')        
    end

end