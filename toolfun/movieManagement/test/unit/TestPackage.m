classdef TestPackage < TestCase
    
    properties
        movie
        package
    end
    
    methods
        function self = TestPackage(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.movie = MovieData();
            self.package = MockPackage(self.movie);
            self.movie.addPackage(self.package);
        end
        
        function tearDown(self)
            delete(self.movie);
            delete(self.package);
        end
        
        %% Tests
        function testGetPackage(self)
            assertEqual(self.movie.getPackage(1), self.package);
        end
        
        
        function testGetSinglePackageIndex(self)
            iPack = self.movie.getPackageIndex(self.package);
            assertEqual(iPack, 1);
        end
        
        function testGetMultiplePackageIndex(self)
            package2 = MockPackage(self.movie);
            self.movie.addPackage(package2);
            iProc = self.movie.getPackageIndex(package2, Inf);
            assertEqual(iProc, [1 2]);
        end
        
        function testDeletePackage(self)
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
        end
        
        function testDeleteInvalidPackage(self)
            % Delete process object
            delete(self.package);
            assertFalse(self.movie.getPackage(1).isvalid);
            
            % Delete process using deletePackage method
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
        end

        %% Tests
        function testCreateDefaultProcess(self)
            self.package.createDefaultProcess(1);
            assertEqual(self.movie.getPackage(1).getProcess(1),...
                self.movie.getProcess(1));
            assertTrue(isa(self.package.getProcess(1),...
                self.package.getProcessClassNames{1}));

        end
    end
end
