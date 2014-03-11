classdef TestProcess < handle
    properties
        movie
        process
    end
    
    methods
        %% Set up and tear down methods
        function setUp(self)
            self.process = MockProcess(self.movie);
            self.movie.addProcess(self.process);
        end
        
        function tearDown(self)
            delete(self.movie);
            delete(self.process);
        end
        
        function package = setUpPackage(self)
            package = MockPackage(self.movie);
            self.movie.addPackage(package);
            package.setProcess(1, self.process);
        end
        
        % Basic tests
        function testIsProcess(self)
            assertTrue(Process.isProcess(class(self.movie.getProcess(1))));
            assertFalse(Process.isProcess('InvalidProcessClass'));
        end
        
        %% Get functions
        function testGetProcess(self)
            assertEqual(self.movie.getProcess(1), self.process);
        end
        
        
        function testGetSingleProcessIndex(self)
            iProc = self.movie.getProcessIndex(self.process);
            assertEqual(iProc, 1);
        end
        
        function testGetMultipleProcessIndex(self)
            process2 = MockProcess(self.movie);
            self.movie.addProcess(process2);
            iProc = self.movie.getProcessIndex(process2, Inf);
            assertEqual(iProc, [1 2]);
        end
        
        %% Delete function
        function testDeleteProcessByIndex(self)
            % Delete process and test deletion by index
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
        end
        
        function testDeleteProcessByProcess(self)
            % Delete process and test deletion by process
            self.movie.deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
        end
        
        function testDeleteSameClassProcessByObject(self)
            % Duplicate process and test deletion by index
            process2 = MockProcess(self.movie);
            self.movie.addProcess(process2);
            
            self.movie.deleteProcess(1);
            assertEqual(self.movie.processes_, {process2});
        end
        
        function testDeletePackageLinkedProcessByIndex(self)
            % Link process to package and test deletion by index
            
            self.setUpPackage();
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
            assertTrue(isempty(self.movie.getPackage(1).getProcess(1)));
        end
        
        function testDeletePackageLinkedProcessByObject(self)
            % Link process to package and test deletion by index
            
            self.setUpPackage();
            self.movie.deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
            assertTrue(isempty(self.movie.getPackage(1).getProcess(1)));
        end
        
        function testDeleteMultiPackageLinkedProcessByIndex(self)
            % Link process to package and test deletion by index
            
            self.setUpPackage();
            self.setUpPackage();
            self.movie.deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
            assertTrue(isempty(self.movie.getPackage(1).getProcess(1)));
            assertTrue(isempty(self.movie.getPackage(2).getProcess(1)));
        end
        
        function testDeleteMultiPackageLinkedProcessByObject(self)
            % Link process to package and test deletion by index
            
            self.setUpPackage();
            self.setUpPackage();
            self.movie.deleteProcess(self.process);
            assertTrue(isempty(self.movie.processes_));
            assertTrue(isempty(self.movie.getPackage(1).getProcess(1)));
            assertTrue(isempty(self.movie.getPackage(2).getProcess(1)));
        end
        
        function testDeleteUnlinkedProcess(self)
            % Delete process and test deletion
            f= @() self.movie.deleteProcess(MockProcess(self.movie));
            assertExceptionThrown(f ,'');
        end
        
        function testDeleteInvalidProcessByIndex(self)
            % Delete process object
            delete(self.process);
            assertFalse(self.movie.getProcess(1).isvalid);
            
            % Delete process using deleteProcess method
            self.movie.deleteProcess(1);
            assertTrue(isempty(self.movie.processes_));
        end
        
        function testReplaceProcess(self)
            process2 = MockProcess(self.movie);
            self.movie.replaceProcess(1, process2);
            
            % Replace process
            assertEqual(self.movie.getProcess(1), process2);
            assertFalse(self.process.isvalid);
        end
        
        function testReplaceProcess2(self)
            process2 = MockProcess(self.movie);
            self.movie.replaceProcess(self.process, process2);
            
            % Replace process
            assertEqual(self.movie.getProcess(1), process2);
            assertFalse(self.process.isvalid);
        end
    end
end
