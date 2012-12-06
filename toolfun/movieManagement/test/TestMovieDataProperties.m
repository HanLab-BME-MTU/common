classdef TestMovieDataProperties < TestCase
    
    properties
        movie
        validProperties =  {'timeInterval_','numAperture_',...
            'magnification_','camBitdepth_'};
        validValues =  {1,1.4,100,14};;
        
    end
    
    methods
        function self = TestMovieDataProperties(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.movie= MovieData();
        end
        
        function tearDown(self)
            delete(self.movie);
        end
        
        %% Tests
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
    end
end
