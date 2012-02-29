classdef TestMovieList < TestCase
    
    properties
        movieList
        movieListPath = fullfile(getenv('HOME'),'MovieTest');
    end
    
    methods
        function self = TestMovieList(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.movieList=TestHelperMovieObject.setUpMovieList(self.movieListPath);
        end
        
        
        function tearDown(self)
            delete(self.movieList);
            rmdir(self.movieListPath,'s');
        end
        
        %% Tests
%         function testGetProcessIndex(self)
%             TestHelperMovieObject.testGetProcessIndex(self.movieList);
%         end
%         
%         
%         function testProcessCreation(self)
%             TestHelperMovieObject.testProcessCreation(self.movieList);
%         end
%         
%         function testPackageCreation(self)
%             TestHelperMovieObject.testPackageCreation(self.movieList);
%         end
        
        
        function testRelocate(self)
            relocatedMovieListPath = TestHelperMovieObject.relocateMovie(self.movieList);
            
            % Load the relocated movie
            relocatedMovieList=MovieData.load(fullfile(relocatedMovieListPath,self.movieList.getFilename),false);
            
            % Test movie paths
            assertEqual(relocatedMovieList.outputDirectory_,relocatedMovieListPath);
            assertEqual(relocatedMovieList.getPath,relocatedMovieListPath);
            
            % Test channel paths
            assertEqual(relocatedMovieList.getMovies{1}.getPath,relocatedMovieListPath);
            
            % Test process/packages relocation            
            rmdir(relocatedMovieListPath,'s');
        end
        
        function testLoad(self)
            newMovieList=MovieList.load(self.movieList.getFullPath);
            assertEqual(self.movieList,newMovieList);
        end
        
        function testSanityCheck(self)
            assertFalse(isempty(self.movieList.getMovies));
            assertTrue(isa(self.movieList.getMovies{1},'MovieData'));
        end
        
        
        function testClass(self)
            assertTrue(isa(self.movieList,'MovieList'));
        end
    end
end
