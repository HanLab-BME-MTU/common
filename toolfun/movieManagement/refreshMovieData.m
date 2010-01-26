function movieData = refreshMovieData(movieData)

%movieData = refreshMovieData(movieData or movieArray)
%
%Loads the most recent version of the input movieData or movie array
%
%Hunter Elliott, 3/2009

wasSingle = false;

%Convert single movie data so it's processed the same
if ~iscell(movieData)
    movieData = {movieData};
    wasSingle = true;
end

nMovies = length(movieData);

for j = 1:nMovies

    %Make sure the OS is right before attempting to load
    movieData{j} = convertMovieDataOS(movieData{j});

    tmp = load([movieData{j}.analysisDirectory filesep 'movieData.mat'],'movieData');


    movieData{j} = tmp.movieData;

    %Convert again incase the loaded moviedata was from linux
    movieData{j} = convertMovieDataOS(movieData{j});

end

%Convert back to structure if single movieData was input
if wasSingle
    movieData = movieData{1};
end