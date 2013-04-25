classdef TestHelperMovieObject < handle
    methods (Static)
        %% Set up and tear down methods
        function movieList = setUpMovieList(path,movie)
            if nargin<2, movie=TestHelperMovieObject.setUpMovie(path);end
            if ~exist(path,'dir'), mkdir(path); end
            movieList = MovieList(movie,path);
            movieList.setPath(path);
            movieList.setFilename('movieList.mat');
            movieList.sanityCheck;
        end
        
        function movie = setUpMovie(path,channel)
            if nargin<2, channel=TestHelperMovieObject.setUpChannel(path);end
            if ~exist(path,'dir'), mkdir(path); end
            movie = MovieData(channel,path);
            movie.setPath(path);
            movie.setFilename('movieData.mat');
        end
        
        function channel = setUpChannel(path,imSize,nFrames, format)
            if nargin<4, format = 'double'; end
            if nargin<3, nFrames=1;end
            if nargin<2,imSize=[100 200]; end
            if ~exist(path,'dir'), mkdir(path); end
            for i=1:nFrames
                imwrite(zeros(imSize, format),fullfile(path,['test_' num2str(i) '.tif']));
            end
            channel=Channel(path);
        end
        
        function relocatedMoviePath= relocateMovie(movieObject)
            % Copy movie in new location
            relocatedMoviePath = [movieObject.getPath '_relocated'];
            copyfile(movieObject.getPath,relocatedMoviePath);
            rmdir(movieObject.getPath, 's')
        end        
    end
end