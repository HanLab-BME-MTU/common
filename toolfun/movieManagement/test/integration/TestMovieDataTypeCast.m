classdef TestMovieDataTypeCast < TestCase

    properties
        movie
        moviePath = fullfile(getenv('HOME'),'MovieTest');
        movieName = 'movieData.mat';
        imSize = [100 200];
        nFrames = 1;            
    end

    methods
        function self = TestMovieDataTypeCast(name)
            self = self@TestCase(name);
        end

        %% Set up and tear down methods
        function setUp(self)
            channel = TestHelperMovieObject.setUpChannel(self.moviePath,self.imSize,self.nFrames);
            self.movie = TestHelperMovieObject.setUpMovie(self.moviePath,channel);
        end
        
        function tearDown(self)
            delete(self.movie);
            rmdir(self.moviePath, 's');
        end
        
        %% Tests    
        function testUINT8(self)
            TestHelperMovieObject.setUpChannel(...
                self.movie.channels_(1).channelPath_,...
                self.imSize,self.nFrames, 'uint8');
            I = self.movie.getChannel(1).loadImage(1);
            assertEqual(class(I), 'uint8');
        end

        function testUINT16(self)
            TestHelperMovieObject.setUpChannel(...
                self.movie.channels_(1).channelPath_,...
                self.imSize,self.nFrames, 'uint16');
            I = self.movie.getChannel(1).loadImage(1);
            assertEqual(class(I), 'uint16');
        end
    end
end
