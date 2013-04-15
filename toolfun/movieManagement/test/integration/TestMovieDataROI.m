classdef TestMovieDataROI < TestCase
    
    properties
        movie
        moviePath = fullfile(getenv('HOME'),'MovieTest');
        movieName = 'movieData.mat';
        roiFolder = 'ROI';
        roiName = 'roiMovie.mat';
        roiMaskName = 'roiMask.tif';
        imSize = [100 200];
        nFrames = 1;
    end
    
    methods
        function self = TestMovieDataROI(name)
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
        
        % ROI tests
        function roiFullPaths = setupROI(self, nROIs)
            
            roiFullPaths = cell(nROIs, 1);
            % Create ROI folder
            for i = 1 : nROIs
                roiPath = fullfile(self.moviePath, self.roiFolder, num2str(i));
                mkdir(roiPath);
                
                % Create ROI mask
                roiMask = true(self.movie.imSize_);
                roiMaskFullPath = fullfile(roiPath, self.roiMaskName);
                imwrite(roiMask, roiMaskFullPath);
                
                % Create and save ROI
                self.movie.addROI(roiMaskFullPath, roiPath);
                self.movie.rois_(i).setPath(roiPath);
                self.movie.rois_(i).setFilename(self.roiName);
                self.movie.rois_(i).sanityCheck;
                roiFullPaths{i} = self.movie.rois_(i).getFullPath();
            end
        end
        
        function testAddROI(self)
            % Create ROI movie
            roiFullPaths = self.setupROI(1);
            
            % Reload ROI Movie and test components
            roiMovie = MovieData.load(roiFullPaths{1});
            assertEqual(roiMovie.parent_, self.movie);
            assertEqual(roiMovie.getROIMask(), true(self.movie.imSize_));
            assertEqual(roiMovie.getAncestor(), self.movie);
            assertEqual(roiMovie.getAncestor().getDescendants(), roiMovie);
        end
        
        function testAddMultipleROIs(self)
            % Create 3 ROIs
            nROIs = 5;
            roiFullPaths = self.setupROI(nROIs);
            assertEqual(numel(self.movie.rois_), nROIs);
            
            % Reload ROI Movie and test components
            roiMovies(1, nROIs) = MovieData();
            for i = 1 : nROIs
                roiMovies(i) = MovieData.load(roiFullPaths{i});
                assertEqual(roiMovies(i).getAncestor(), self.movie);
            end
            
            % Test getAncestor/getDescendants methods
            assertEqual(self.movie.getDescendants(), roiMovies);
        end
        
        function testDeleteROI(self)
            % Create ROI movie
            roiFullPaths = self.setupROI(2);
            assertEqual(numel(self.movie.rois_),2);
            
            % Delete create ROI
            self.movie.deleteROI(1);
            assertEqual(numel(self.movie.rois_),1);
            self.movie.save;
            
            % Test ROI has been deleted
            roiMovie1 = load(roiFullPaths{1});
            assertFalse(roiMovie1.MD.isvalid)
            
            % Test ROI has been deleted
            roiMovie2 = MovieData.load(roiFullPaths{2});
            assertEqual(roiMovie2.getAncestor(), self.movie);
            
            % Test main movie has been saved without ROI
            reloadedMovie = MovieData.load(self.movie.getFullPath);
            assertEqual(reloadedMovie.getDescendants(), roiMovie2)
            assertEqual(numel(reloadedMovie.rois_),1)
        end
        
        
        function testRelocateROI(self)
            % Add ROI & save
            self.setupROI(1)
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
        
        function testDeleteSharedProcess(self)
            
            % Create process in main movie
            self.movie.addProcess(ThresholdProcess(self.movie));
            
            % Create ROI movies
            self.setupROI(2);
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
        
        function testDeleteSharedPackage(self)
            
            % Create process in main movie
            self.movie.addPackage(SegmentationPackage(self.movie));
            
            % Create ROI movies
            self.setupROI(2);
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
        
        function testReplaceSharedProcess(self)
            % Create process
            self.movie.addProcess(ThresholdProcess(self.movie));
            oldprocess = self.movie.getProcess(1);
            
            % Create ROI movies
            self.setupROI(2);
            
            % Replace process
            self.movie.replaceProcess(1, MaskRefinementProcess(self.movie))
            assertTrue(isa(self.movie.getProcess(1),'MaskRefinementProcess'));
            assertTrue(isa(self.movie.rois_(1).getProcess(1),'MaskRefinementProcess'));
            assertTrue(isa(self.movie.rois_(2).getProcess(1),'MaskRefinementProcess'));
            assertFalse(oldprocess.isvalid);
        end
    end
end
