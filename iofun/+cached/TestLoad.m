classdef TestLoad < TestCase
    properties
        matfile
        ascii_file
    end
    methods
        function self = TestLoad(name)
            self = self@TestCase(name);
        end
        function setUp(self)
            self.matfile = [ tempname '.mat' ];
            S.background = [];
            S.psfSigma = 1.0793;
            S.movieInfo(200) = struct('xCoord',5,'yCoord',10,'amp',5);
            save(self.matfile,'-struct','S');
            
            self.ascii_file = [ tempname '.txt' ];
            p = rand(10);
            q = rand(10);
            save(self.ascii_file,'p','q','-ascii');
        end
        function tearDown(self)
            cached.load('-clear');
            delete(self.matfile);
            delete(self.ascii_file);
        end
        function testMatFileName(self,varargin)
            S = cached.load(self.matfile,varargin{:});
            fileStruct = load(self.matfile,varargin{:});
            assertEqual(fileStruct,S);
            
            S = cached.load(self.matfile,'-mat',varargin{:});
            assertEqual(fileStruct,S);
        end
        function testAsciiFileName(self,varargin)
            S = cached.load(self.ascii_file,varargin{:});
            fileStruct = load(self.ascii_file,varargin{:});
            assertEqual(fileStruct,S);
            
            S = cached.load(self.ascii_file,'-ascii',varargin{:});
            assertEqual(fileStruct,S);
        end
        function testMatVariables(self)
            self.testMatFileName('background');
            self.testMatFileName('psfSigma','movieInfo');
            self.testMatFileName;
        end
        function testReset(self,varargin)
            [S, wasCached] = cached.load(self.matfile,varargin{:});
            assert(~wasCached);
            [S, wasCached] = cached.load(self.matfile,varargin{:});
            assert(wasCached);
            [S, wasCached] = cached.load(self.matfile,varargin{:},'-reset');
            assert(~wasCached);
            
            S.newVar = struct('a',5,'b',10);
            save(self.matfile,'-struct','S');
            [S, wasCached] = cached.load(self.matfile,varargin{:});
            assert(~isfield(S,'newVar'));
            [S, wasCached] = cached.load(self.matfile,varargin{:},'-reset','newVar');
            assert(isfield(S,'newVar'));
        end
        function testUseCache(self,varargin)
            [~, wasCached] = cached.load(self.matfile,varargin{:},'-useCache',true);
            assert(~wasCached);
            [~, wasCached] = cached.load(self.matfile,varargin{:},'-useCache',true);
            assert(wasCached);
            [S, wasCached] = cached.load(self.matfile,varargin{:},'-useCache',false);
            assert(~wasCached);
            
            S.newVar = struct('a',5,'b',10);
            save(self.matfile,'-struct','S');
            [S, wasCached] = cached.load(self.matfile,varargin{:},'-useCache',true);
            assert(~isfield(S,'newVar'));
            [S, wasCached] = cached.load(self.matfile,varargin{:},'-useCache',false,'newVar');
            assert(isfield(S,'newVar'));
        end
        function testTerminalUseCache(self)
            self.testUseCache('psfSigma','movieInfo');
        end
        function testTerminalReset(self)
            self.testReset('background');
        end
        function testClear(self)
            [~, wasCached] = cached.load(self.matfile);
            assert(~wasCached);
            [~, wasCached] = cached.load(self.ascii_file);
            assert(~wasCached);
            
            [~, wasCached] = cached.load(self.matfile);
            assert(wasCached);
            [~, wasCached] = cached.load(self.ascii_file);
            assert(wasCached);
            
            cached.load('-clear');
            
            [~, wasCached] = cached.load(self.matfile);
            assert(~wasCached);
            [~, wasCached] = cached.load(self.ascii_file);
            assert(~wasCached);
        end
        function testClearException(self)
            try
                [~, wasCached] = cached.load(self.matfile,'background','-clear');
                % an exception should have been thrown since -clear must be
                % the only argument
                assert(false);
            catch
                % exception thrown
            end
        end
    end
end