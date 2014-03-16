classdef TestMovieDataProperties < TestCase & TestLibrary
    
    properties
        timeInterval =  1
        numAperture = 1.4
        magnification = 100
        camBitdepth = 14
    end
    
    methods
        function self = TestMovieDataProperties(name)
            self = self@TestCase(name);
        end
        
        function setUp(self)
            self.setUpMovieData();
        end
        
        function tearDown(self)
            tearDown@TestLibrary(self);
        end
        
        %% Individual property tests
        function testSetValidTimeInterval(self)
            self.movie.timeInterval_ = self.timeInterval;
            assertEqual(self.movie.timeInterval_, self.timeInterval);
        end
        
        function testSetInvalidTimeInterval(self)
            f= @() set(self.movie, 'timeInterval_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidNumericalAperture(self)
            self.movie.numAperture_ = self.numAperture;
            assertEqual(self.movie.numAperture_, self.numAperture);
        end
        
        function testSetInvalidNumericalAperture(self)
            f= @() set(self.movie, 'numAperture_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidMagnification(self)
            self.movie.magnification_ = self.magnification;
            assertEqual(self.movie.magnification_, self.magnification);
        end
        
        function testSetInvalidMagnification(self)
            f= @() set(self.movie, 'magnification_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        function testSetValidCamBitdepth(self)
            self.movie.camBitdepth_ = self.camBitdepth;
            assertEqual(self.movie.camBitdepth_, self.camBitdepth);
        end
        
        function testSetInvalidCamBitdepth(self)
            f= @() set(self.movie, 'camBitdepth_', 0);
            assertExceptionThrown(f,'lccb:set:invalid');
        end
        
        %% Multi
        function testSetMultipleProperties(self)
            properties = {'timeInterval_', 'numAperture_',...
                'magnification_', 'camBitdepth_'};
            values = {self.timeInterval, self.numAperture,...
                self.magnification, self.camBitdepth};
            set(self.movie, properties, values);
            for i = 1 : numel(properties)
                assertEqual(self.movie.(properties{i}), values{i});
            end
        end
        
        function testMultiSetProperties(self)
            properties = {'timeInterval_', 'numAperture_',...
                'magnification_', 'camBitdepth_'};
            values = {self.timeInterval, self.numAperture,...
                self.magnification,self.camBitdepth};
            set(self.movie, properties, values);
            set(self.movie, properties, values);
            for i = 1 : numel(properties)
                assertEqual(self.movie.(properties{i}), values{i});
            end
        end
    end
end
