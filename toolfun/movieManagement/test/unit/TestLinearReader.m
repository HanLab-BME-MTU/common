classdef TestLinearReader < TestProxyReader
    methods
        function self = TestLinearReader(name,reader)
            if(nargin < 2)
                reader = MockReader;
            end
            self = self@TestProxyReader(name,LinearReader(reader));
        end
        function checkFcnLinear(proxyFcn,varargin)
               rSize = [ proxy.reader.getSizeC , proxy.reader.getSizeT, proxy.reader.getSizeZ ];
               linFcn = @(p,c,t,z) proxyFcn(p,sub2ind(rSize,c,t,z));
               checkFcnToZ(linFcn,varargin{:});
        end

    end
end
