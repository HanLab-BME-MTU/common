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
            movie=TestHelperMovieObject.setUpMovie(self.movieListPath);
            movie.sanityCheck();
            self.movieList=TestHelperMovieObject.setUpMovieList(self.movieListPath, movie);
        end
        
        
        function tearDown(self)
            delete(self.movieList);
            if isdir(self.movieListPath)
                rmdir(self.movieListPath,'s');
            end
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
        
        %% Process/package deletion
        
        function testDeleteSingleProcess(self)            
            % Create process
            self.movieList.addProcess(SignalPreprocessingProcess(self.movieList,'',struct()));
            assertEqual(numel(self.movieList.processes_),1);
            
            % Create process and test deletion
            self.movieList.deleteProcess(1);
            assertTrue(isempty(self.movieList.processes_));
        end
        
        function testReplaceSingleProcess(self)
            % Create process
            self.movieList.addProcess(SignalPreprocessingProcess(self.movieList,'',struct()));
            oldprocess = self.movieList.getProcess(1);
            
            % Replace process
            self.movieList.replaceProcess(1, SignalProcessingProcess(self.movieList,'',struct()))
            assertTrue(isa(self.movieList.getProcess(1), 'SignalProcessingProcess'));
            assertFalse(oldprocess.isvalid);
        end
        
        function testDeleteSinglePackage(self)
            % Create package
            self.movieList.addPackage(IntegratorPackage(self.movieList));
            assertEqual(numel(self.movieList.packages_),1);
            
            % Delete package and test deletion
            self.movieList.deletePackage(1);
            assertTrue(isempty(self.movieList.packages_));
        end
    end
end
