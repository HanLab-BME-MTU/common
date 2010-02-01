function status = checkMovieMaskMovie(movieData)

% status = checkMovieMaskMovie(movieData)
% 
% This function checks whether a mask movie has been successfully created
% using makeMaskMovie.m and returns status = true if so, and status = false
% if not.
% 
% Hunter Elliott, 10/2009
% 


status = false;

if isfield(movieData,'masks') && isfield(movieData.masks,'movie') ...
        && isfield(movieData.masks.movie,'status') && movieData.masks.status == 1 ...
        && isfield(movieData.masks.movie,'fileName') && ...
        exist([movieData.analysisDirectory filesep movieData.masks.movie.fileName],'file')
    status = true;
    
end