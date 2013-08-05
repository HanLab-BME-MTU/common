classdef TestBFMovieData < TestMovieData & TestCase
    
    properties
        fakename = 'test.fake';
    end
    
    methods
        function self = TestBFMovieData(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.setUp@TestMovieData();
        end
        
        function tearDown(self)
            self.tearDown@TestMovieData();
        end
        
        function setUpMovie(self)
            filename = fullfile(self.path, self.fakename);
            fid = fopen(filename, 'w');
            fclose(fid);
            
            self.movie = MovieData.load(filename);
        end
        
        function checkChannelPaths(self)
            for i = 1 : self.nChan
                assertEqual(self.movie.getChannel(i).channelPath_,...
                    fullfile(self.path, self.fakename))
            end
        end
        
        %% Typecasting tests
        function checkPixelType(self, classname)
            if strcmp(classname, 'single'),
                self.fakename = 'test&pixelType=float.fake';
            else
                self.fakename = ['test&pixelType=' classname '.fake'];
            end
            self.setUpMovie();
            I = self.movie.getChannel(1).loadImage(1);
            assertTrue(isa(I, classname));
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
        
        %% Dimensions tests
        function testSizeXY(self)
            self.fakename = 'test&sizeX=256&sizeY=256.fake';
            self.imSize = [256 256];
            self.setUpMovie()
            self.checkMovie();
        end
        
        function testSizeC(self)
            self.fakename = 'test&sizeC=4.fake';
            self.nChan = 4;
            self.setUpMovie()
            self.checkMovie();
        end
        
        function testSizeT(self)
            self.fakename = 'test&sizeT=256.fake';
            self.nFrames = 256;
            self.setUpMovie()
            self.checkMovie();
        end
    end
end
