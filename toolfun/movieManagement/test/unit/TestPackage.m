classdef TestPackage < handle
    
    properties
        movie
        package
    end
    
    methods
        %% Set up and tear down methods
        function setUp(self)
            self.package = MockPackage(self.movie);
            self.movie.addPackage(self.package);
        end
        
        function process = setUpProcess(self)
            process = MockProcess(self.movie);
            self.movie.addProcess(process);
            self.package.setProcess(1, process);
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
        
        %% deletePackage tests
        function testDeletePackageByIndex(self)
            % Delete package by index
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
        end
        
        function testDeletePackageByObject(self)
            % Test package deletion by object
            self.movie.deletePackage(self.package);
            assertTrue(isempty(self.movie.packages_));
        end
        
        function testDeleteSameClassPackageByIndex(self)
            % Duplicate package class and test deletion by index
            package2 = MockPackage(self.movie);
            self.movie.addPackage(package2);
            
            self.movie.deletePackage(1);
            assertEqual(self.movie.packages_, {package2});
        end
        
        function testDeleteSameClassPackageByObject(self)
            % Duplicate package class and test deletion by object
            package2 = MockPackage(self.movie);
            self.movie.addPackage(package2);
            
            self.movie.deletePackage(self.package);
            assertEqual(self.movie.packages_, {package2});
        end
        
        function testDeleteLoadedPackageByIndex(self)
            % Link process to package and test deletion by index
            
            process = self.setUpProcess();
            self.movie.deletePackage(1);
            assertEqual(self.movie.processes_, {process});
        end
        
        function testDeleteLoadedPackageByObject(self)
            % Link process to package and test deletion by object
            
            process = self.setUpProcess();
            self.movie.deletePackage(self.package);
            assertEqual(self.movie.processes_, {process});
        end
        
        function testDeleteInvalidPackage(self)
            % Delete process object
            delete(self.package);
            assertFalse(self.movie.getPackage(1).isvalid);
            
            % Delete process using deletePackage method
            self.movie.deletePackage(1);
            assertTrue(isempty(self.movie.packages_));
        end
        
        %% createDefaultProcess tests
        function testCreateDefaultProcess(self)
            self.package.createDefaultProcess(1);
            assertEqual(self.movie.getPackage(1).getProcess(1),...
                self.movie.getProcess(1));
            assertTrue(isa(self.package.getProcess(1),...
                self.package.getProcessClassNames{1}));
            
        end
    end
end
