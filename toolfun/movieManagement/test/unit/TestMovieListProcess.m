classdef TestMovieListProcess < TestCase & TestProcess
    
    methods
        function self = TestMovieListProcess(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.movie = MovieList();
            setUp@TestProcess(self);
        end
        
        function tearDown(self)
            tearDown@TestProcess(self);
        end
    end
end
