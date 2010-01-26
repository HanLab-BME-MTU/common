function movieData = revertMovieData(movieData)

%This loads the previous saved movieData, reverting the moviedata back to
%it's last state. 
%
%Hunter Elliott, 3/2009

%Check if it is a single movie data or an array
if ~iscell(movieData)
    wasSingle = true;
    movieData = {movieData};
else
    wasSingle = false;
end

nMovies = length(movieData);

for j = 1:nMovies


    %movieData{j} = convertMovieDataOS(movieData{j});

    %Load the old moviedata
    tmp = load([movieData{j}.analysisDirectory filesep 'movieData.old'],'-mat');
    movieData{j} = tmp.movieData;

end

%Return as a structure if necessary
if wasSingle
    movieData = movieData{1};
end