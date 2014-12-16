classdef TestMovieData < TestMovieObject
    
    properties
        imSize = [512 512];
        nFrames = 1;
        nChan = 1;
        zSize = 1;
        
        % ROI properties
        roiFolder = 'ROI';
        roiName = 'roiMovie.mat';
        roiMaskName = 'roiMask.tif';
    end
    
    methods
        %% Set up and tear down methods
        function setUp(self)
            self.setUp@TestMovieObject();
        end
        
        function tearDown(self)
            delete(self.movie);
            self.tearDown@TestMovieObject();
        end
        
        %% SanityCheck test
        function checkDimensions(self)
            assertEqual(numel(self.movie.channels_), self.nChan);
            assertEqual(self.movie.imSize_, self.imSize);
            assertEqual(self.movie.nFrames_, self.nFrames);
            assertEqual(self.movie.zSize_, self.zSize);
        end
        
        function checkMovie(self)
            assertTrue(isa(self.movie,'MovieData'));
            self.checkDimensions()
            self.checkChannelPaths();
        end
        
        function rois = setUpROIs(self, nROIs, roiMask)
            
            if nargin < 3, roiMask = true(self.movie.imSize_); end
            % Create ROI folder
            rois(nROIs, 1) = MovieData();
            for i = 1 : nROIs
                roiPath = fullfile(self.movie.getPath(),...
                    [self.roiFolder '_' num2str(i)]);
                mkdir(roiPath);
                
                % Create ROI mask
                roiMaskFullPath = fullfile(roiPath, self.roiMaskName);
                imwrite(roiMask, roiMaskFullPath);
                
                % Create and save ROI
                rois(i) = self.movie.addROI(roiMaskFullPath, roiPath);
                rois(i).setPath(roiPath);
                rois(i).setFilename(self.roiName);
                rois(i).sanityCheck;
            end
        end
        
        %% Tests
        function testSimple(self)
            self.setUpMovie();
            self.checkMovie();
        end
        
        function testRelocate(self)
            self.setUpMovie();
            
            % Perform movie relocation
            moviePath = self.movie.getPath();
            movieName = self.movie.getFilename();
            oldPath = self.path;
            self.relocate();
            
            % Load the relocated movie
            newPath = relocatePath(moviePath, oldPath, self.path);
            newFullPath = fullfile(newPath, movieName);
            self.movie = MovieData.load(newFullPath, false);
            self.checkMovie();
            
            % Test movie paths
            assertEqual(self.movie.outputDirectory_, newPath);
            assertEqual(self.movie.getPath, newPath);
        end
        
        function testLoadAbsolutePath(self)
            % Test MovieData loading from absolute path
            self.setUpMovie();
            self.movie = MovieData.load(self.movie.getFullPath());
            self.checkMovie();
        end
        
        function testLoadRelativePath(self)
            % Test MovieData loading from relative path
            here = pwd;
            self.setUpMovie();
            absolutePath = self.movie.getPath();
            cd(self.movie.getPath());
            relativePath = fullfile('.', self.movie.getFilename());
            self.movie = MovieData.load(relativePath);
            self.checkMovie();
            assertEqual(self.movie.getPath(), absolutePath);
            cd(here);
        end
        
        function testLoadSymlink(self)
            % Test MovieData loading from a symlink
            
            if ispc, return; end
            self.setUpMovie();
            absolutePath = self.movie.getPath();
            symlinkPath = self.createSymlink(absolutePath);
            symlinkFullPath = fullfile(symlinkPath, self.movie.getFilename());
            self.movie = MovieData.load(symlinkFullPath);
            self.checkMovie();
            assertEqual(self.movie.getPath(), absolutePath);            
        end
        
        %% ROI
        function testSimpleROI(self)
            self.setUpMovie();
            rois = self.setUpROIs(1);
            roiMovieFullPath = rois(1).getFullPath();
            
            % Test ROI has been deleted
            self.movie = MovieData.load(roiMovieFullPath);
            self.checkMovie();
        end
        
        function testGetROIMaskDefault(self)
            self.setUpMovie();
            
            % Read the movie and check the roi mask is read and cached
            roiMask = true([self.imSize self.nFrames]);
            assertTrue(isempty(self.movie.roiMask));
            assertEqual(self.movie.getROIMask(), roiMask);
            assertEqual(self.movie.roiMask, roiMask);
            
            % Reload the movie and check the roimask is transient
            self.movie = MovieData.load(self.movie.getFullPath());
            assertTrue(isempty(self.movie.roiMask));
        end
        
        function testGetROIMask(self)
            self.setUpMovie();
            roiMask = true([self.imSize self.nFrames]);
            roiMask(1:end/2, 1:end/2, :) = false;
            roi = self.setUpROIs(1, roiMask);
            
            % Read the movie and check the roi mask is read and cached
            assertTrue(isempty(roi.roiMask));
            assertEqual(roi.getROIMask(), roiMask);
            assertEqual(roi.roiMask, roiMask);
            
            % Reload the movie and check the roimask is transient
            roi = MovieData.load(roi.getFullPath());
            assertTrue(isempty(roi.roiMask));
        end
        
        function testDeleteROI(self)
            % Create ROI movie
            self.setUpMovie();
            rois = self.setUpROIs(2);
            assertEqual(numel(self.movie.rois_), 2);
            roiMovie1FullPath = rois(1).getFullPath();
            roiMovie2FullPath = rois(2).getFullPath();
            
            % Delete create ROI
            self.movie.deleteROI(1);
            assertEqual(numel(self.movie.rois_), 1);
            self.movie.sanityCheck();
            
            % Test ROI has been deleted
            roiMovie1 = load(roiMovie1FullPath);
            assertFalse(roiMovie1.MD.isvalid)
            
            % Test ROI has been deleted
            roiMovie2 = MovieData.load(roiMovie2FullPath);
            assertEqual(numel(roiMovie2.getAncestor().getDescendants), 1);
            
            % Test main movie has been saved without ROI
            self.movie = MovieData.load(self.movie.getFullPath);
            assertEqual(numel(self.movie.getDescendants()), 1)
        end
        
        function testRelocateROI(self)
            % Add ROI & save
            self.setUpMovie();
            self.setUpROIs(1);
            
            self.movie.sanityCheck();
            
            % Perform movie relocation
            moviePath = self.movie.getROI(1).getPath();
            movieName = self.movie.getROI(1).getFilename();
            roiOutputDirectory = self.movie.getROI(1).outputDirectory_;
            roiMaskPath = self.movie.getROI(1).roiMaskPath_;
            oldPath = self.path;
            self.relocate();
            
            % Load the relocated ROI
            newPath = relocatePath(moviePath, oldPath, self.path);
            newFullPath = fullfile(newPath, movieName);
            self.movie = MovieData.load(newFullPath, false);
            
            % Test movie paths
            self.checkMovie();
            assertEqual(self.movie.outputDirectory_,...
                relocatePath(roiOutputDirectory, oldPath, self.path));
            assertEqual(self.movie.roiMaskPath_,...
                relocatePath(roiMaskPath, oldPath, self.path));
        end
        
        function testSharedProcess(self)
            self.setUpMovie();
            self.setUpProcess();
            self.setUpROIs(1);
            self.movie.sanityCheck();
            
            % Test package
            self.movie = MovieData.load(self.movie.getFullPath());
            assertEqual(self.movie.getProcess(1),...
                self.movie.getROI(1).getProcess(1));
        end
        
        function testSharedPackage(self)
            self.setUpMovie();
            self.setUpPackage();
            self.setUpROIs(1);
            self.movie.sanityCheck();
            
            % Test package
            self.movie = MovieData.load(self.movie.getFullPath());
            assertEqual(self.movie.getPackage(1),...
                self.movie.getROI(1).getPackage(1));
        end
    end
    
    methods (Abstract)
        setUpMovie(self)
        checkChannelPaths(self)
    end
end
