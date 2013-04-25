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

    end
end
