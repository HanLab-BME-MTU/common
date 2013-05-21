classdef TestOmeroMovieData < TestCase & TestOmero
    
    methods
        function self = TestOmeroMovieData(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.movie = MovieData();
            setUp@TestOmero(self);
        end
        
        function tearDown(self)
            tearDown@TestOmero(self);
        end
    end
end
