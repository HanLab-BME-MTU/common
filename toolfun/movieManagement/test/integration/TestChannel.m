classdef TestChannel < TestCase
    
    properties
        channel
        channelPath = fullfile(getenv('HOME'),'Desktop','ChannelTest');
        imSize=[100 200];
    end
    
    methods
        function self = TestChannel(name)
            self = self@TestCase(name);
        end
        
        function setUp(self)
            self.channel = TestHelperMovieObject.setUpChannel(self.channelPath,self.imSize,1);
        end
        
        function tearDown(self)
            delete(self.channel);
            rmdir(self.channelPath,'s');
        end
        
        function testSanityCheck(self)
            [width,height,nf]= self.channel.sanityCheck;
            assertEqual(width,self.imSize(2));
            assertEqual(height,self.imSize(1));
            assertEqual(nf,1);
        end
        
        function testPath(self)
            assertEqual(self.channel.channelPath_,self.channelPath);
        end
        
    end
end
