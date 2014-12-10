classdef TestMovieObject < TestLibrary
    
    properties
        path
    end
    
    methods
        %% Set up and tear down methods
        function setUp(self)
            uuid = char(java.util.UUID.randomUUID().toString());
            self.path = fullfile(self.tmpdir, uuid);
            if ~exist(self.path, 'dir'), mkdir(self.path); end
            [~, f] = fileattrib(self.path);
            self.path = f.Name;
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
