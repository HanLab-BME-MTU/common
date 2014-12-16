classdef TestMovieList < TestMovieObject & TestCase
    
    properties
        movieList
        nMovies = 1;
    end
    
    methods
        function self = TestMovieList(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.setUp@TestMovieObject();
        end
        
        function setUpMovieList(self)
            filename = fullfile(self.path, 'test.fake');
            fid = fopen(filename, 'w');
            fclose(fid);
            movie = MovieData.load(filename);
            
            self.movieList = MovieList(movie, self.path);
            self.movieList.setPath(self.path);
            self.movieList.setFilename('movieList.mat');
            self.movieList.sanityCheck();
        end
        
        function tearDown(self)
            delete(self.movieList);
            self.tearDown@TestMovieObject();
        end
        
        %% SanityCheck test
        function checkMovieList(self)
            assertTrue(isa(self.movieList,'MovieList'));
            assertFalse(isempty(self.movieList.getMovies));
            
            for i = 1 : self.nMovies
                assertTrue(isa(self.movieList.getMovie(i), 'MovieData'));
            end
        end
        
        
        %% Tests
        
        function testSimple(self)
            self.setUpMovieList();
            self.checkMovieList();
        end
        
        function testRelocate(self)
            self.setUpMovieList();
            
            % Perform relocation
            movieListPath = self.movieList.getPath();
            movieListName = self.movieList.getFilename();
            oldPath = self.path;
            self.relocate();
            
            % Load the relocated movie
            newPath = relocatePath(movieListPath, oldPath, self.path);
            newFullPath = fullfile(newPath, movieListName);
            self.movieList = MovieList.load(newFullPath, false);
            self.checkMovieList();
            
            % Test movie paths
            assertEqual(self.movieList.outputDirectory_, newPath);
            assertEqual(self.movieList.getPath, newPath);
        end
        
        function testLoadAbsolutePath(self)
            % Test MovieList loading from relative path
            
            self.setUpMovieList();
            self.movieList = MovieList.load(self.movieList.getFullPath());
            self.checkMovieList;
        end
        
        function testLoadRelativePath(self)
            % Test MovieList loading from relative path
            here = pwd;
            self.setUpMovieList();
            absolutePath = self.movieList.getPath();
            cd(self.movieList.getPath());
            relativePath = fullfile('.', self.movieList.getFilename());
            self.movieList = MovieList.load(relativePath);
            self.checkMovieList();
            assertEqual(self.movieList.getPath(), absolutePath);
            cd(here);
        end
        
        function testLoadSymlink(self)
            % Test MovieList loading from a symlink
            
            if ispc, return; end
            self.setUpMovieList();
            absolutePath = self.movieList.getPath();
            symlinkPath = self.createSymlink(absolutePath);
            symlinkFullPath = fullfile(symlinkPath, self.movieList.getFilename());
            self.movieList = MovieList.load(symlinkFullPath);
            self.checkMovieList();
            assertEqual(self.movieList.getPath(), absolutePath);            
        end
    end
end
