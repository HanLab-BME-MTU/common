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
    end
end
