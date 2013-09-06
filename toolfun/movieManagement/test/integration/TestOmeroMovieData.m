classdef TestOmeroMovieData < TestMovieData & TestCase
    
    properties
        client
        session
        id
    end
    
    methods
        function self = TestOmeroMovieData(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.setUp@TestMovieData();
            [self.client, self.session] = loadOmero();
        end
        
        function tearDown(self)
            self.client.closeSession();
            self.tearDown@TestMovieData();
        end
        
        function setUpMovie(self, type)
            pixelsService = self.session.getPixelsService();
            self.id = pixelsService.createImage(self.imSize(2),...
                self.imSize(1), 1, self.nFrames,...
                self.nChan, type, 'test image', 'TestOmeroMovieData');
            
            
            self.movie = MovieData.load(self.session, self.id);
        end
        
        function checkChannels(self)
            for i = 1 : self.nChan
                assertEqual(self.movie.getChannel(i).channelPath_,...
                    fullfile(self.path, self.fakename))
            end
        end
        
        %% Typecasting tests
        function checkPixelType(self, type)
            self.setUpMovie(type);
            self.checkMovie();
        end
        
        function testINT8(self)
            self.checkPixelType('int8');
        end
        
        function testUINT8(self)
            self.checkPixelType('uint8');
        end
        
        function testINT16(self)
            self.checkPixelType('int16');
        end
        
        function testUINT16(self)
            self.checkPixelType('uint16');
        end
        
        function testUINT32(self)
            self.checkPixelType('uint32');
        end
        
        function testSINGLE(self)
            self.checkPixelType('single');
        end
        
        function testDOUBLE(self)
            self.checkPixelType('double');
        end        
    end
end
