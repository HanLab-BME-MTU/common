classdef TestBFMovieData < TestCase

    properties
        movie
        moviePath = fullfile(getenv('HOME'),'MovieTest');
        movieName = 'movieData.mat';
        imSize = [10 50];
        nFrames = 3;
        nChan = 2;
    end

    methods
        function self = TestBFMovieData(name)
            self = self@TestCase(name);
        end

        %% Set up and tear down methods
        function setUp(self)
        end
        
        function tearDown(self)
            delete(self.movie);
            rmdir(self.moviePath,'s');
        end
        
        %% SanityCheck test            
        function testSanityCheck(self)
            self.movie=TestHelperMovieObject.setUpBFMovie(self.moviePath,...
                self.imSize, self.nChan, self.nFrames);
            assertEqual(self.movie.channels_(1).owner_,self.movie);
            assertEqual(self.movie.imSize_,self.imSize);
            assertEqual(self.movie.nFrames_,self.nFrames);
            assertEqual(numel(self.movie.channels_),self.nChan);
        end

        %% Typecasting tests
        function testUINT8(self)
            self.movie=TestHelperMovieObject.setUpBFMovie(self.moviePath,...
                self.imSize, self.nChan, self.nFrames, 'uint8');
            I = self.movie.channels_(1).loadImage(1);            
            assertEqual(I, zeros(self.imSize, 'uint8'));
        end
        
        function testUINT16(self)
            self.movie=TestHelperMovieObject.setUpBFMovie(self.moviePath,...
                self.imSize, self.nChan, self.nFrames, 'uint16');
            I = self.movie.channels_(1).loadImage(1);            
            assertEqual(I, zeros(self.imSize, 'uint16'));
        end
        function testUINT32(self)
            self.movie=TestHelperMovieObject.setUpBFMovie(self.moviePath,...
                self.imSize, self.nChan, self.nFrames, 'uint32');
            I = self.movie.channels_(1).loadImage(1);            
            assertEqual(I, zeros(self.imSize, 'uint32'));
        end
        function testSINGLE(self)
            self.movie=TestHelperMovieObject.setUpBFMovie(self.moviePath,...
                self.imSize, self.nChan, self.nFrames, 'single');
            I = self.movie.channels_(1).loadImage(1);            
            assertEqual(I, zeros(self.imSize, 'uint32'));
        end
        function testUINT64(self)
            self.movie=TestHelperMovieObject.setUpBFMovie(self.moviePath,...
                self.imSize, self.nChan, self.nFrames, 'double');
            I = self.movie.channels_(1).loadImage(1);            
            assertEqual(I, zeros(self.imSize, 'uint64'));
        end
    end
end
