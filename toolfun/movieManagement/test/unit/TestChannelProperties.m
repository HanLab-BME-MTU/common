classdef TestChannelProperties < TestCase
    
    properties
        channel
        emissionWavelength = 500;
        excitationWavelength = 600
        imageTypes = Channel.getImagingModes()
        exposureTime = 100
        fluorophores = Channel.getFluorophores()
    end
    
    methods
        function self = TestChannelProperties(name)
            self = self@TestCase(name);
        end
        
        function setUp(self)
            self.channel = Channel();
        end
        
        function tearDown(self)
            delete(self.channel);
        end
        
        %% Class test
        function testClass(self)
            assertTrue(isa(self.channel,'Channel'));
        end
        
        %% Individual property tests
        function testSetValidEmissionWavelength(self)
            self.channel.emissionWavelength_ = self.emissionWavelength;
            assertEqual(self.channel.emissionWavelength_, self.emissionWavelength);
        end
        
        function testSetInvalidEmissionWavelength(self)
            f= @() set(self.channel, 'emissionWavelength_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidExcitationWavelength(self)
            self.channel.excitationWavelength_ = self.excitationWavelength;
            assertEqual(self.channel.excitationWavelength_, self.excitationWavelength);
        end
        
        function testSetInvalidExcitationWavelength(self)
            f= @() set(self.channel, 'excitationWavelength_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidExposureTime(self)
            self.channel.exposureTime_ = self.exposureTime;
            assertEqual(self.channel.exposureTime_, self.exposureTime);
        end
        
        function testSetInvalidExposureTime(self)
            f= @() set(self.channel, 'exposureTime_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        function testSetValidiImageType(self)
            for i = 1 : numel(self.imageTypes)
                self.channel = Channel();
                self.channel.imageType_ = self.imageTypes{i};
                assertEqual(self.channel.imageType_, self.imageTypes{i});
            end
        end
        
        function testSetInvalidImageType(self)
            f= @() set(self.channel, 'imageType_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidFlurophores(self)
            for i = 1 : numel(self.fluorophores)
                self.channel = Channel();
                self.channel.fluorophore_ = self.fluorophores{i};
                assertEqual(self.channel.fluorophore_, self.fluorophores{i});
            end
        end
        
        function testSetInvalidFlurophore(self)
            f= @() set(self.channel, 'fluorophore_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
    end
end
