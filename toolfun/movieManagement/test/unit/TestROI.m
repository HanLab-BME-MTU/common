classdef TestROI < TestCase
    
    properties
        movie
        roi = MovieData.empty(0,1);
        nRois = 5
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
        
        %% Single ROI tests
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
    end
end
