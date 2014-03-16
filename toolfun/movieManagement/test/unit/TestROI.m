classdef TestROI < TestCase & TestLibrary
        
    methods
        function self = TestROI(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.setUpMovieData();
        end
        
        function tearDown(self)
            tearDown@TestLibrary(self);
        end       
        
        %% ROI methods tests
        function testAddROI(self)
            self.rois = self.movie.addROI('','');
            assertEqual(self.movie.rois_, self.rois);
            assertEqual(self.rois.parent_, self.movie);
        end
        
        function testGetROI(self)
            rois = self.setUpRois();
            for i = 1 : self.nRois
                assertEqual(self.movie.getROI(i), rois(i));
            end
        end
        
        function testGetAncestor(self)
            self.setUpRois();
            for i = 1 : self.nRois
                assertEqual(self.rois(i).getAncestor(), self.movie);
            end
        end
        
        function testGetDescendants(self)
            self.setUpRois();
            assertEqual(self.movie.getDescendants(), self.rois);
        end
        
        function testDeleteROI(self)
            rois = self.setUpRois();
            self.movie.deleteROI(1, false);
            assertEqual(self.movie.rois_, rois(2 : end));
            assertFalse(rois(1).isvalid);
        end
        
        function testDeleteROIs(self)
            self.setUpRois();
            self.movie.deleteROI(1 : self.nRois - 1, false);
            assertEqual(self.movie.rois_, self.rois(self.nRois));
            for i = 1 : self.nRois -1
                assertFalse(self.rois(i).isvalid);
            end
        end
        
        % Shared process/package tests
        function testSharedChannels(self)
            channels = [Channel() Channel()];
            self.movie = MovieData(channels, '');
            self.setUpRois();
            for i = 1 : self.nRois
                assertEqual(self.movie.getROI(i).channels_, self.movie.channels_);
            end
        end
        
        function testSharedMetadata(self)
            self.movie.pixelSize_ = 100;
            self.setUpRois();
            for i = 1 : self.nRois
                assertEqual(self.movie.getROI(i).channels_, self.movie.channels_);
            end
        end
        
        function testSharedProcess(self)
            self.setUpProcess();
            self.setUpRois();
            assertEqual(self.movie.processes_, self.processes);
            for i = 1: self.nRois
                assertEqual(self.movie.getROI(i).processes_, self.processes);
            end
        end
        
        function testUnsharedProcess(self)
            self.setUpRois();
            self.setUpProcess();
            assertEqual(self.movie.processes_, self.processes);
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessByIndex(self)
            self.setUpProcess();
            self.setUpRois();
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessByObject(self)
            process = self.setUpProcess();
            self.setUpRois();
            self.movie.deleteProcess(process);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessFromROIByIndex(self)
            process = self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(self.nRois).deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessFromROIByObject(self)
            process = self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(self.nRois).deleteProcess(process);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testUnlinkSharedProcessFromROIByIndex(self)
            process = self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(1).unlinkProcess(1);
            assertEqual(self.movie.processes_, {process});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).processes_, {process});
            end
            assertTrue(isempty(self.movie.getROI(1).processes_));
        end
        
        function testUnlinkSharedProcessFromROIByObject(self)
            process = self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(1).unlinkProcess(process);
            assertEqual(self.movie.processes_, {process});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).processes_, {process});
            end
            assertTrue(isempty(self.movie.getROI(1).processes_));
        end
        
        function testReplaceSharedProcess(self)
            self.setUpProcess();
            self.setUpRois();
            newprocess = MockProcess(self.movie);
            self.movie.replaceProcess(1, newprocess);
            assertEqual(self.movie.getProcess(1), newprocess);
            for i = 1: self.nRois
                assertEqual(self.movie.getROI(i).getProcess(1), newprocess);
            end
        end
        
        function testSharedPackage(self)
            package = self.setUpPackage();
            self.setUpRois();
            assertEqual(self.movie.getPackage(1), package);
            for i = 1: self.nRois
                assertEqual(self.movie.getROI(i).packages_, {package});
            end
        end
        
        function testUnsharedPackage(self)
            self.setUpRois();
            package = self.setUpPackage();
            assertEqual(self.movie.getPackage(1), package);
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testDeleteSharedPackageByIndex(self)
            self.setUpPackage();
            self.setUpRois();
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testDeleteSharedPackageByObject(self)
            package = self.setUpPackage();
            self.setUpRois();
            self.movie.deletePackage(package);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testDeleteSharedPackageFromROIByIndex(self)
            self.setUpPackage();
            self.setUpRois();
            self.movie.getROI(self.nRois).deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testDeleteSharedPackageFromROIByObject(self)
            package = self.setUpPackage();
            self.setUpRois();
            self.movie.getROI(self.nRois).deletePackage(package);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testUnlinkSharedPackageFromROIByIndex(self)
            package = self.setUpPackage();
            self.setUpRois();
            self.movie.getROI(1).unlinkPackage(1);
            assertEqual(self.movie.packages_, {package});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).packages_, {package});
            end
            assertTrue(isempty(self.movie.getROI(1).packages_));
        end
        
        function testUnlinkSharedPackageFromROIByObject(self)
            package = self.setUpPackage();
            self.setUpRois();
            self.movie.getROI(1).unlinkPackage(package);
            assertEqual(self.movie.packages_, {package});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).packages_, {package});
            end
            assertTrue(isempty(self.movie.getROI(1).packages_));
        end
        
        %% cleanupROIPackages integration tests
        function testCleanupROIPackagesNoKeep(self)
            self.setUpPackage(true);
            self.setUpRois();
            cleanupROIPackages(self.movie, 'MockPackage');
            for i = 1: numel(self.movie.rois_)
                assertTrue(isempty(self.movie.getROI(i).packages_));
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testCleanupROIPackagesKeep(self)
            
            [package, process] = self.setUpPackage(true);
            self.setUpRois();
            cleanupROIPackages(self.movie, 'MockPackage', 1);
            
            % Tests
            for i = 1: numel(self.movie.rois_)
                roiPackage = self.movie.getROI(i).getPackage(1);
                assertTrue(isa(roiPackage, 'MockPackage'));
                assertFalse(isequal(package, roiPackage));
                assertEqual(roiPackage.owner_, self.movie.getROI(i));
                assertEqual(roiPackage.getProcess(1), process);
            end
        end
        
        function testCleanupROIPackagesFromChild(self)
            
            self.setUpPackage(true);
            self.setUpRois();
            cleanupROIPackages(self.movie.getROI(1), 'MockPackage');
            for i = 1: numel(self.movie.rois_)
                assertTrue(isempty(self.movie.getROI(i).packages_));
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testCleanupROIPackagesNoROI(self)
            
            package = self.setUpPackage(true);
            cleanupROIPackages(self.movie, 'MockPackage');
            assertEqual(self.movie.packages_, {package})
        end
    end
end
