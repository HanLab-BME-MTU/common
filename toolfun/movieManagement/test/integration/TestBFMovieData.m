classdef TestBFMovieData < TestCase
    
    properties
        movie
        path = fullfile(getenv('HOME'),'MovieTest');
        name = 'movieData.mat';
        imSize = [10 50];
        nFrames = 3;
        nChan = 2;
    end
    
    methods
        function self = TestBFMovieData(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function movie = setUpBFMovie(self, format)
            if nargin < 2, format = 'double'; end
            if ~exist(self.path,'dir'), mkdir(self.path); end
            bfsave(zeros(self.imSize(1), self.imSize(2), self.nChan, self.nFrames, format),...
                fullfile(self.path, 'test.ome.tiff'), 'XYCTZ');
            movie = bfImport(fullfile(self.path, 'test.ome.tiff'));
        end
        
        function tearDown(self)
            delete(self.movie);
            rmdir(self.path,'s');
        end
        
        %% SanityCheck test
        function testSanityCheck(self)
            self.movie= self.setUpBFMovie();
            assertEqual(self.movie.channels_(1).owner_, self.movie);
            assertEqual(self.movie.imSize_, self.imSize);
            assertEqual(self.movie.nFrames_, self.nFrames);
            assertEqual(numel(self.movie.channels_), self.nChan);
        end
        
        %% Typecasting tests
        function checkPixelType(self, classname)
            self.movie = self.setUpBFMovie(classname);
            I = self.movie.getChannel(1).loadImage(1);
            assertEqual(I, zeros(self.imSize, classname));
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
