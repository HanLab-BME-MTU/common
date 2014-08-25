classdef TestProxyReader < TestCase
    properties
        proxy
        reader
    end
    methods
        function self = TestProxyReader(name,proxy,reader)
            self = self@TestCase(name);
            if(nargin < 2)
                self.proxy = ProxyReader(MockReader);
            else
                self.proxy = proxy;
            end
            if(nargin < 3)
                self.reader = self.proxy.reader;
            end
        end
        function checkFcn(self,fcnProxy,fcnReader)
            if(nargin < 3)
                fcnReader = fcnProxy;
            end
            for c = 1 : self.reader.getSizeC
                assertEqual( fcnProxy(self.proxy,c) , ... 
                             fcnReader(self.reader,c) ) ;
            end
        end
        function checkFcnToZ(self,fcnProxy,fcnReader)
            if(nargin < 3)
                fcnReader = fcnProxy;
            end
            for c = 1 : self.reader.getSizeC
                for t = 1 : self.reader.getSizeT(c)
                    for z = 1 : self.reader.getSizeZ(c)
                        assertEqual( fcnProxy(self.proxy,c,t,z) , ... 
                                     fcnReader(self.reader,c,t,z) );
                    end
                end
            end
        end
        function testGetSizeX(self)
            self.checkFcn(@getSizeX);
        end
        function testGetSizeY(self)
            self.checkFcn(@getSizeY);
        end
        function testGetSizeZ(self)
            self.checkFcn(@getSizeZ);
        end
        function testGetSizeT(self)
            self.checkFcn(@getSizeT);
        end
        function testGetSizeC(self)
            assertEqual(self.proxy.getSizeC,self.reader.getSizeC);
        end
        function testGetBitDepth(self)
            self.checkFcn(@getBitDepth);
        end
        function testGetImageFileNames(self)
            self.checkFcn(@getImageFileNames);
        end
        function testGetChannelNames(self)
            self.checkFcn(@getChannelNames);
        end
        function testLoadImage(self)
            self.checkFcnToZ(@loadImage);
        end
        function testLoadStack(self)
            self.checkFcnToZ(@loadStack);
        end
    end
end
