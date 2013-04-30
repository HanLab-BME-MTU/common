classdef TestMovieDataPackage < TestCase & TestPackage
    
    methods
        function self = TestMovieDataPackage(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.movie = MovieData();
            setUp@TestPackage(self);
        end
        
        function tearDown(self)
            tearDown@TestPackage(self);
        end
    end
end
