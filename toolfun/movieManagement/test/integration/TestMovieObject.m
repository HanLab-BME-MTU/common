classdef TestMovieObject < handle
    
    properties
        path = fullfile(getenv('HOME'), 'MovieTest');
    end
    
    methods
        %% Set up and tear down methods
        function setUp(self)
            if ~exist(self.path, 'dir'), mkdir(self.path); end
        end
        
        function tearDown(self)
            rmdir(self.path, 's');
        end
        
        %% Library methods
        function relocate(self)
            relocatedPath = [self.path '_relocated'];
            movefile(self.path, relocatedPath);
            self.path = relocatedPath;
        end
    end
end
