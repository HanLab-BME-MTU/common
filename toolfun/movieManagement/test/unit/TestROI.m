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
            for i = 1 : self.nRois
                self.roi(i) = self.movie.addROI('','');
            end
        end
        
        function tearDown(self)
            delete(self.movie);
        end
        
        %% ROI methods tests
        function testGetROI(self)
            for i = 1 : self.nRois
                assertEqual(self.movie.getROI(i), self.roi(i));
            end
        end
        
        function testGetAncestor(self)
            for i = 1 : self.nRois
                assertEqual(self.roi(i).getAncestor(), self.movie);
            end
        end
        
        function testGetDescendants(self)
            assertEqual(self.movie.getDescendants(), self.roi);
        end
        
        function testDeleteROI(self)
            self.movie.deleteROI(1, false);
            assertEqual(self.movie.rois_, self.roi(2 : self.nRois));
            assertFalse(self.roi(1).isvalid);
        end
        
        function testDeleteROIs(self)
            self.movie.deleteROI(1 : self.nRois - 1, false);
            assertEqual(self.movie.rois_, self.roi(self.nRois));
            for i = 1 : self.nRois -1
                assertFalse(self.roi(i).isvalid);
            end
        end
        
        % Shared process/package tests
        function setUpSharedProcess(self)
            self.process = MockProcess(self.movie);
            self.movie.addProcess(self.process);
            
            self.nRois = self.nRois + 1;
            self.roi(self.nRois) = self.movie.addROI('','');
        end
        
        function setUpSharedPackage(self, loaded)
            self.package = MockPackage(self.movie);
            self.movie.addPackage(self.package);
            
            if nargin > 1 && loaded
                self.package.createDefaultProcess(1);
                self.process = self.package.getProcess(1);
            end
            
            self.nRois = self.nRois + 1;
            self.roi(self.nRois) = self.movie.addROI('','');
        end
        
        function testSharedProcess(self)
            self.setUpSharedProcess();
            assertEqual(self.movie.getProcess(1), self.process);
            for i = 1: self.nRois - 1
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
            assertEqual(self.movie.getROI(self.nRois).getProcess(1), self.process);
        end
        
        function testDeleteSharedProcessByIndex(self)
            self.setUpSharedProcess();
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessByObject(self)
            self.setUpSharedProcess();
            self.movie.deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessFromROIByIndex(self)
            self.setUpSharedProcess();
            self.movie.getROI(self.nRois).deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testDeleteSharedProcessFromROIByObject(self)
            self.setUpSharedProcess();
            self.movie.getROI(self.nRois).deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testUnlinkSharedProcessFromROIByIndex(self)
            self.setUpSharedProcess();
            self.movie.getROI(self.nRois).unlinkProcess(1);
            assertEqual(self.movie.processes_, {self.process});
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testUnlinkSharedProcessFromROIByObject(self)
            self.setUpSharedProcess();
            self.movie.getROI(self.nRois).unlinkProcess(self.process);
            assertEqual(self.movie.processes_, {self.process});
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testReplaceSharedProcess(self)
            self.setUpSharedProcess();
            newprocess = MockProcess(self.movie);
            self.movie.replaceProcess(1, newprocess);
            assertEqual(self.movie.getProcess(1), newprocess);
            for i = 1: self.nRois - 1
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
            assertEqual(self.movie.getROI(self.nRois).getProcess(1), newprocess);
        end
        
        function testSharedPackage(self)
            self.setUpSharedPackage();
            assertEqual(self.movie.getPackage(1), self.package);
            for i = 1: self.nRois - 1
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
            assertEqual(self.movie.getROI(self.nRois).getPackage(1), self.package);
        end
        
        function testDeleteSharedPackageByIndex(self)
            self.setUpSharedPackage();
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testDeleteSharedPackageByObject(self)
            self.setUpSharedPackage();
            self.movie.deletePackage(self.package);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testDeleteSharedPackageFromROIByIndex(self)
            self.setUpSharedPackage();
            self.movie.getROI(self.nRois).deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testDeleteSharedPackageFromROIByObject(self)
            self.setUpSharedPackage();
            self.movie.getROI(self.nRois).deletePackage(self.package);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testUnlinkSharedPackageFromROIByIndex(self)
            self.setUpSharedPackage();
            self.movie.getROI(self.nRois).unlinkPackage(1);
            assertEqual(self.movie.packages_, {self.package});
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        function testUnlinkSharedPackageFromROIByObject(self)
            self.setUpSharedPackage();
            self.movie.getROI(self.nRois).unlinkPackage(self.package);
            assertEqual(self.movie.packages_, {self.package});
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
        
        
        %% cleanupROIPackages integration tests
        function testCleanupROIPackagesNoKeep(self)
            
            self.setUpSharedPackage(true);
            cleanupROIPackages(self.movie, 'MockPackage');
            for i = 1: numel(self.movie.rois_)
                assertTrue(isempty(self.movie.getROI(i).packages_));
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
        
        function testCleanupROIPackagesKeep(self)
            
            self.setUpSharedPackage(true);
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
            
            self.setUpSharedPackage(true);
            cleanupROIPackages(self.movie.getROI(1), 'MockPackage');
            for i = 1: numel(self.movie.rois_)
                assertTrue(isempty(self.movie.getROI(i).packages_));
                assertTrue(isempty(self.movie.getROI(i).processes_));
            end
        end
    end
end
