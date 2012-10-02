classdef TestMovieData < TestCase

    properties
        movie
        moviePath = fullfile(getenv('HOME'),'MovieTest');
        movieName = 'movieData.mat';
        roiFolder = 'ROI';
        roiName = 'roiMovie.mat';
        roiMaskName = 'roiMask.tif';
        imSize = [100 200];
        nFrames = 1;
        validProperties =  {'timeInterval_','numAperture_',...
                'magnification_','camBitdepth_'};
        validValues =  {1,1.4,100,14};;
            
    end

    methods
        function self = TestMovieData(name)
            self = self@TestCase(name);
        end

        %% Set up and tear down methods
        function setUp(self)
            channel = TestHelperMovieObject.setUpChannel(self.moviePath,self.imSize,self.nFrames);
            self.movie=TestHelperMovieObject.setUpMovie(self.moviePath,channel);
        end
        
        function tearDown(self)
            delete(self.movie);
            rmdir(self.moviePath,'s');
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
        
        function testSetInvalidProperties(self)
            TestHelperMovieObject.testSetInvalidProperties(self.movie,self.validProperties);
        end
        
        function testMultiSetProperties(self)
            TestHelperMovieObject.testMultiSetProperties(self.movie,self.validProperties,self.validValues);
        end
        
        function testSetMultipleProperties(self)
            TestHelperMovieObject.testSetMultipleProperties(self.movie,self.validProperties,self.validValues);
        end
        
        
        function testSetIndividualProperties(self)
            TestHelperMovieObject.testSetIndividualProperties(self.movie,self.validProperties,self.validValues);
        end
        
        function testSanityCheck(self)
            assertEqual(self.movie.channels_.owner_,self.movie);
            assertEqual(self.movie.imSize_,self.imSize);
            assertEqual(self.movie.nFrames_,1);
        end
          
        function testClass(self)
            assertTrue(isa(self.movie,'MovieData'));
        end
        
        % ROI tests
        function testAddROI(self)
            % Create ROI folder
            roiPath = fullfile(self.moviePath, self.roiFolder);
            mkdir(roiPath);

            % Create ROI mask
            roiMask = true(self.movie.imSize_); 
            roiMaskFullPath = fullfile(roiPath, self.roiMaskName);
            imwrite(roiMask, roiMaskFullPath);
            
            % Create and save ROI
            self.movie.addROI(roiMaskFullPath, roiPath);
            self.movie.rois_(end).setPath(roiPath);
            self.movie.rois_(end).setFilename(self.roiName);
            self.movie.rois_(end).sanityCheck;
            
            % Reload ROI Movie and test components
            roiMovie = MovieData.load(fullfile(roiPath, self.roiName));
            assertEqual(roiMovie.parent_, self.movie);
            assertEqual(roiMovie.getROIMask(), true(self.movie.imSize_));
        end
        
        % ROIS test
    end
end
