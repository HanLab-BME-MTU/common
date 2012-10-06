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
        function roiFullPath = setupROI(self)
            % Create ROI folder
            nROIs = numel(self.movie.rois_);
            roiPath = fullfile(self.moviePath, self.roiFolder, num2str(nROIs+1));
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
            roiFullPath = self.movie.rois_(end).getFullPath();
        end
        
        function testAddROI(self)
            % Create ROI movie
            roiFullPath = self.setupROI();
            
            % Reload ROI Movie and test components
            roiMovie = MovieData.load(roiFullPath);
            assertEqual(roiMovie.parent_, self.movie);
            assertEqual(roiMovie.getROIMask(), true(self.movie.imSize_));
            assertEqual(roiMovie.getAncestor(), self.movie);
            assertEqual(roiMovie.getAncestor().getDescendants(), roiMovie);
        end
        
        function testAddMultipleROIs(self)
            % Create 3 ROIs
            roiFullPath{1} = self.setupROI();
            roiFullPath{2} = self.setupROI();
            roiFullPath{3} = self.setupROI();
            assertEqual(numel(self.movie.rois_), 3);
            
            % Reload ROI Movie and test components
            roiMovie1 = MovieData.load(roiFullPath{1});
            roiMovie2 = MovieData.load(roiFullPath{2});
            roiMovie3 = MovieData.load(roiFullPath{3});
            
            % Test getAncestor/getDescendants methods
            assertEqual(roiMovie1.getAncestor(), self.movie);
            assertEqual(roiMovie2.getAncestor(), self.movie);
            assertEqual(roiMovie3.getAncestor(), self.movie);
            assertEqual(self.movie.getDescendants(),...
                [roiMovie1, roiMovie2, roiMovie3]);

        end
        
        function testDeleteROI(self)
            % Create ROI movie
            roiFullPath{1} = self.setupROI();
            roiFullPath{2} = self.setupROI();
            assertEqual(numel(self.movie.rois_),2);
            
            % Delete create ROI
            self.movie.deleteROI(1);            
            assertEqual(numel(self.movie.rois_),1);
            self.movie.save;
            
            % Test ROI has been deleted
            roiMovie1 = load(roiFullPath{1});
            assertFalse(roiMovie1.MD.isvalid)
            
            % Test ROI has been deleted
            roiMovie2 = MovieData.load(roiFullPath{2});
            assertEqual(roiMovie2.getAncestor(), self.movie);

            % Test main movie has been saved without ROI
            reloadedMovie = MovieData.load(self.movie.getFullPath);
            assertEqual(reloadedMovie.getDescendants(), roiMovie2)
            assertEqual(numel(reloadedMovie.rois_),1)
        end
        
                
        function testRelocateROI(self)
            % Add ROI & save
            self.setupROI()
            self.movie.save;
            roiOutputDirectory = self.movie.rois_(1).outputDirectory_;
            roiMaskPath = self.movie.rois_(1).roiMaskPath_;

             % Load the relocated movie
            relocatedMoviePath = TestHelperMovieObject.relocateMovie(self.movie);
            relocatedMovie=MovieData.load(fullfile(relocatedMoviePath,self.movie.getFilename),false);
            
            % Test movie paths
            relocatedROI = relocatedMovie.rois_(1);
            relocatedROIOutputDir = relocatePath(roiOutputDirectory,...
                self.movie.getPath, relocatedMoviePath);
            relocatedROIMaskPath = relocatePath(roiMaskPath,...
                self.movie.getPath, relocatedMoviePath);
            assertEqual(relocatedROI.outputDirectory_, relocatedROIOutputDir);
            assertEqual(relocatedROI.roiMaskPath_, relocatedROIMaskPath);

            % Remove relocated movie
            rmdir(relocatedMoviePath,'s');
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
        
        
        function testDeleteSharedProcess(self)
            
            % Create process in main movie
            self.movie.addProcess(ThresholdProcess(self.movie));
            
            % Create ROI movies
            self.setupROI();
            self.setupROI();
            assertEqual(numel(self.movie.processes_),1);
            assertEqual(numel(self.movie.rois_(1).processes_),1);
            assertEqual(numel(self.movie.rois_(2).processes_),1);
            
            % Delete process from ROI 1
            self.movie.rois_(1).deleteProcess(1);  
            assertEqual(numel(self.movie.processes_),0);
            assertEqual(numel(self.movie.rois_(1).processes_),0);
            assertEqual(numel(self.movie.rois_(2).processes_),0);
            
            % Reload movie and check process has been deleted in the parent
            % and in the ROI movies
            self.movie.save();
            reloadedMovie = MovieData.load(self.movie.getFullPath);
            assertEqual(numel(reloadedMovie.processes_),0);
            assertEqual(numel(reloadedMovie.rois_(1).processes_),0);
            assertEqual(numel(reloadedMovie.rois_(2).processes_),0);

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

        
        function testDeleteSharedPackage(self)
            
            % Create process in main movie
            self.movie.addPackage(SegmentationPackage(self.movie));
            
            % Create ROI movies
            self.setupROI();
            self.setupROI();
            assertEqual(numel(self.movie.packages_),1);
            assertEqual(numel(self.movie.rois_(1).packages_),1);
            assertEqual(numel(self.movie.rois_(2).packages_),1);
            
            % Delete process
            self.movie.rois_(1).deletePackage(1);  
            assertEqual(numel(self.movie.packages_),0);
            assertEqual(numel(self.movie.rois_(1).packages_),0);
            assertEqual(numel(self.movie.rois_(2).packages_),0);
            
            % Reload movie and check process has been deleted in the parent
            % and in the ROI movies
            self.movie.save();
            reloadedMovie = MovieData.load(self.movie.getFullPath);
            assertEqual(numel(reloadedMovie.packages_),0);
            assertEqual(numel(reloadedMovie.rois_(1).packages_),0);
            assertEqual(numel(reloadedMovie.rois_(2).packages_),0);

        end
        
    end
end
