classdef TestChannel < TestCase

    properties
        channel
        channelPath = fullfile(getenv('HOME'),'Desktop','ChannelTest');
        imSize=[100 200];
        validProperties={'emissionWavelength_','excitationWavelength_','imageType_',...
            'exposureTime_','fluorophore_'};
        validValues={500,600,'TIRF',100,'cfp'};
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
        
        function testSetInvalidProperties(self)
            TestHelperMovieObject.testSetInvalidProperties(self.channel,self.validProperties);
        end
        
        function testMultiSetProperties(self)
            TestHelperMovieObject.testMultiSetProperties(self.channel,self.validProperties,self.validValues);
        end
        
        function testSetMultipleProperties(self)
            TestHelperMovieObject.testSetMultipleProperties(self.channel,self.validProperties,self.validValues);
        end
        
        
        function testSetIndividualProperties(self)
            TestHelperMovieObject.testSetIndividualProperties(self.channel,self.validProperties,self.validValues);
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
        
        
        function testClass(self)
            assertEqual(class(self.channel),'Channel');
        end
          
    end
end
