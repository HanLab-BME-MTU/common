classdef TestMovieDataROIMask < TestCase & TestLibrary
    
    
    methods
        function self = TestMovieDataROIMask(name)
            self = self@TestCase(name);
        end
        
        function setUp(self)
            self.setUpMovieData();
        end
        
        function tearDown(self)
            tearDown@TestLibrary(self);
        end
        
        function testSetInvalidROIMaskPath(self)
            f= @() self.movie.setROIMaskPath(true(2,2));
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidROIMaskPath(self)
            self.movie.setROIMaskPath('/path/to/mask.tif');
            assertEqual(self.movie.roiMaskPath_, '/path/to/mask.tif');
        end
        
        function testSetInvalidROIOmeroId(self)
            f= @() self.movie.setROIOmeroId('/path/to/mask.tif');
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidROIOmeroId(self)
            self.movie.setROIOmeroId(1);
            assertEqual(self.movie.roiOmeroId_, 1);
        end
        

    end
end
