classdef TestROI < TestCase
    
    properties
        movie
        roi = MovieData.empty(0,1);
        nRois = 5
        process
        package
    end
    
    methods
        function self = TestROI(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.movie = MovieData();
        end
        
        function tearDown(self)
            delete(self.movie);
        end
        
        function setUpRois(self)
            for i = 1 : self.nRois
                self.roi(end+1) = self.movie.addROI('','');
            end
        end
        
        function setUpProcess(self)
            self.process = MockProcess(self.movie);
            self.movie.addProcess(self.process);
        end
        
        function setUpPackage(self, loaded)
            self.package = MockPackage(self.movie);
            self.movie.addPackage(self.package);
            
            if nargin > 1 && loaded
                self.package.createDefaultProcess(1);
                self.process = self.package.getProcess(1);
            end
        end
        
        %% ROI methods tests
        function testAddROI(self)
            self.roi = self.movie.addROI('','');
            assertEqual(self.movie.rois_, self.roi);
            assertEqual(self.roi.parent_, self.movie);
        end
        
        function testGetROI(self)
            self.setUpRois();
            for i = 1 : self.nRois
                assertEqual(self.movie.getROI(i), self.roi(i));
            end
        end
        
        function testGetAncestor(self)
            self.setUpRois();
            for i = 1 : self.nRois
                assertEqual(self.roi(i).getAncestor(), self.movie);
            end
        end
        
        function testGetDescendants(self)
            self.setUpRois();
            assertEqual(self.movie.getDescendants(), self.roi);
        end
        
        function testDeleteROI(self)
            self.setUpRois();
            self.movie.deleteROI(1, false);
            assertEqual(self.movie.rois_, self.roi(2 : self.nRois));
            assertFalse(self.roi(1).isvalid);
        end
        
        function testDeleteROIs(self)
            self.setUpRois();
            self.movie.deleteROI(1 : self.nRois - 1, false);
            assertEqual(self.movie.rois_, self.roi(self.nRois));
            for i = 1 : self.nRois -1
                assertFalse(self.roi(i).isvalid);
            end
        end
        
        % Shared process/package tests
        function testSharedProcess(self)
            self.setUpProcess();
            self.setUpRois();
            assertEqual(self.movie.getProcess(1), self.process);
            for i = 1: self.nRois
                assertEqual(self.movie.getROI(i).processes_, {self.process});
            end
        end
        
        function testUnsharedProcess(self)
            self.setUpRois();
            self.setUpProcess();
            assertEqual(self.movie.getProcess(1), self.process);
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
            self.setUpProcess();
            self.setUpRois();
            self.movie.deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessFromROIByIndex(self)
            self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(self.nRois).deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessFromROIByObject(self)
            self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(self.nRois).deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testUnlinkSharedProcessFromROIByIndex(self)
            self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(1).unlinkProcess(1);
            assertEqual(self.movie.processes_, {self.process});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).processes_, {self.process});
            end
            assertTrue(isempty(self.movie.getROI(1).processes_));
        end
        
        function testUnlinkSharedProcessFromROIByObject(self)
            self.setUpProcess();
            self.setUpRois();
            self.movie.getROI(1).unlinkProcess(self.process);
            assertEqual(self.movie.processes_, {self.process});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).processes_, {self.process});
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
            self.setUpPackage();
            self.setUpRois();
            assertEqual(self.movie.getPackage(1), self.package);
            for i = 1: self.nRois
                assertEqual(self.movie.getROI(i).packages_, {self.package});
            end
        end
        
        function testUnsharedPackage(self)
            self.setUpRois();
            self.setUpPackage();
            assertEqual(self.movie.getPackage(1), self.package);
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
            self.setUpPackage();
            self.setUpRois();
            self.movie.deletePackage(self.package);
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
            self.setUpPackage();
            self.setUpRois();
            self.movie.getROI(self.nRois).deletePackage(self.package);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testUnlinkSharedPackageFromROIByIndex(self)
            self.setUpPackage();
            self.setUpRois();
            self.movie.getROI(1).unlinkPackage(1);
            assertEqual(self.movie.packages_, {self.package});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).packages_, {self.package});
            end
            assertTrue(isempty(self.movie.getROI(1).packages_));
        end
        
        function testUnlinkSharedPackageFromROIByObject(self)
            self.setUpPackage();
            self.setUpRois();
            self.movie.getROI(1).unlinkPackage(self.package);
            assertEqual(self.movie.packages_, {self.package});
            for i = 2: self.nRois
                assertEqual(self.movie.getROI(i).packages_, {self.package});
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
            
            self.setUpPackage(true);
            self.setUpRois();
            cleanupROIPackages(self.movie, 'MockPackage', 1);
            
            % Tests
            for i = 1: numel(self.movie.rois_)
                roiPackage = self.movie.getROI(i).getPackage(1);
                assertTrue(isa(roiPackage, 'MockPackage'));
                assertFalse(isequal(self.package, roiPackage));
                assertEqual(roiPackage.owner_, self.movie.getROI(i));
                assertEqual(roiPackage.getProcess(1), self.process);
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
            
            self.setUpPackage(true);
            cleanupROIPackages(self.movie, 'MockPackage');
            assertEqual(self.movie.packages_, {self.package})
        end
    end
end
