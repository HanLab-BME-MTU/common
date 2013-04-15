classdef TestMovieData < TestCase
    
    properties
        movie
        moviePath = fullfile(getenv('HOME'),'MovieTest');
        movieName = 'movieData.mat';
        imSize = [100 200];
        nFrames = 1;
    end
    
    methods
        function self = TestMovieData(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            channel = TestHelperMovieObject.setUpChannel(self.moviePath,self.imSize,self.nFrames);
            self.movie=TestHelperMovieObject.setUpMovie(self.moviePath,channel);
            self.movie.sanityCheck();
        end
        
        function tearDown(self)
            delete(self.movie);
            if isdir(self.moviePath)
                rmdir(self.moviePath,'s');
            end
        end
        
        %% Tests
        function testGetProcessIndex(self)
            TestHelperMovieObject.testGetProcessIndex(self.movie);
        end
        
        function testRelocate(self)
            % Add process + package to test analysis relocation
            self.movie.addProcess(ThresholdProcess(self.movie));
            self.movie.addPackage(SegmentationPackage(self.movie));
            self.movie.sanityCheck;
            self.movie.save;
            
            % Load the relocated movie
            relocatedMoviePath = TestHelperMovieObject.relocateMovie(self.movie);
            relocatedMovie=MovieData.load(fullfile(relocatedMoviePath,self.movie.getFilename),false);
            
            % Test movie paths
            assertEqual(relocatedMovie.outputDirectory_,relocatedMoviePath);
            assertEqual(relocatedMovie.getPath,relocatedMoviePath);
            
            % Test process/packages paths
            assertEqual(fileparts(relocatedMovie.processes_{1}.funParams_.OutputDirectory),relocatedMoviePath);
            assertEqual(fileparts(relocatedMovie.packages_{1}.outputDirectory_),relocatedMoviePath);
            
            rmdir(relocatedMoviePath,'s');
        end
        
        function testLoad(self)
            newProcess= ThresholdProcess(self.movie);
            newPackage= SegmentationPackage(self.movie);
            self.movie.addProcess(newProcess);
            self.movie.addPackage(newPackage);
            self.movie.sanityCheck;
            self.movie.save;
            newmovie=MovieData.load(self.movie.getFullPath);
            assertEqual(self.movie,newmovie);
            assertEqual(newmovie.processes_,{newProcess});
            assertEqual(newmovie.packages_,{newPackage});
        end
        
        function testClass(self)
            assertTrue(isa(self.movie,'MovieData'));
        end
        
        %% Process/package deletion
        function testDeleteSingleProcess(self)
            % Create process
            self.movie.addProcess(ThresholdProcess(self.movie));
            assertEqual(numel(self.movie.processes_),1);
            
            % Delete process and test deletion
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
        end
        
        function testDeleteInvalidProcess(self)
            
            % Create process
            self.movie.addProcess(ThresholdProcess(self.movie));
            assertEqual(numel(self.movie.processes_),1);
            
            % Delete process object
            delete(self.movie.processes_{1});
            assertEqual(numel(self.movie.processes_),1);
            
            % Delete process using deleteProcess method
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
        end
        
        function testDeleteSinglePackage(self)
            % Create package
            self.movie.addPackage(SegmentationPackage(self.movie));
            assertEqual(numel(self.movie.packages_),1);
            
            % Delete package and test deletion
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
        end
        
        function testDeleteInvalidPackage(self)
            
            % Create process in main movie
            self.movie.addPackage(SegmentationPackage(self.movie));
            assertEqual(numel(self.movie.packages_),1);
            
            % Delete package object
            delete(self.movie.packages_{1});
            assertEqual(numel(self.movie.packages_),1);
            
            % Delete package using deletePackage method
            self.movie.deletePackage(1);
            assertEqual(numel(self.movie.packages_),0);
        end
        
        function testReplaceSingleProcess(self)
            % Create process
            self.movie.addProcess(ThresholdProcess(self.movie));
            oldprocess = self.movie.getProcess(1);
            
            % Replace process
            self.movie.replaceProcess(1, MaskRefinementProcess(self.movie))
            assertTrue(isa(self.movie.getProcess(1),'MaskRefinementProcess'));
            assertFalse(oldprocess.isvalid);
        end
    end
end
