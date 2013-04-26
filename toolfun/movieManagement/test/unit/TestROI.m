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
        
        function setUpSharedPackage(self)
            self.package = MockPackage(self.movie);
            self.movie.addPackage(self.package);
            
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
        
        function testDeleteSharedProcess(self)
            self.setUpSharedProcess();
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
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
        
        function testDeleteSharedPackage(self)
            self.setUpSharedPackage();
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
            for i = 1: self.nRois
                assertTrue(isempty(self.movie.getROI(i).packages_));
            end
        end
    end
end
