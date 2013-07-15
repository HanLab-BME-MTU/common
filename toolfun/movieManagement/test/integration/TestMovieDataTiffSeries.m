classdef TestMovieDataTiffSeries < TestMovieData & TestCase
    
    
    methods
        function self = TestMovieDataTiffSeries(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.setUp@TestMovieData();
        end
        
        function tearDown(self)
            self.tearDown@TestMovieData();
        end
        
        function setUpMovie(self, format)
            if nargin < 2,  format = 'uint8'; end
            channel = TestHelperMovieObject.setUpChannel(self.path,self.imSize,...
                self.nFrames, format);
            self.movie=TestHelperMovieObject.setUpMovie(self.path,channel);
            self.movie.sanityCheck();
        end
        
        function checkChannels(self)
            for i = 1 : self.nChan
                assertEqual(self.movie.getChannel(i).channelPath_,self.path);
            end
        end
        
        %% Pixel Type tests
        function testUINT8(self)
            self.setUpMovie('uint8')
            I = self.movie.getChannel(1).loadImage(1);
            assertEqual(class(I), 'uint8');
        end
        
        function testUINT16(self)
            self.setUpMovie('uint16');
            self.movie = MovieData.load(self.movie.getFullPath());
            I = self.movie.getChannel(1).loadImage(1);
            assertEqual(class(I), 'uint16');
        end
        %% Invalid sizeT tests
        function testInvalidNumberFrames(self)
            self.setUpMovie();
            fullPath = self.movie.getFullPath();
            TestHelperMovieObject.setUpChannel(self.movie.getChannel(1).channelPath_,...
                self.imSize(end:-1:1),self.nFrames+1);
            assertExceptionThrown(@() MovieData.load(fullPath), 'MovieData:sanityCheck:nFrames');
        end
        
        function testInvalidImSize(self)
            self.setUpMovie();
            fullPath = self.movie.getFullPath();
            TestHelperMovieObject.setUpChannel(self.movie.getChannel(1).channelPath_,...
                self.imSize/2,self.nFrames);
            assertExceptionThrown(@() MovieData.load(fullPath), 'MovieData:sanityCheck:imSize');
        end
        
    end
end
