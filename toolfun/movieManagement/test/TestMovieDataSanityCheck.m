classdef TestMovieDataSanityCheck < TestCase

    properties
        movie
        moviePath = fullfile(getenv('HOME'),'MovieTest');
        movieName = 'movieData.mat';
        imSize = [100 200];
        nFrames = 3;
        nChan = 2;
    end

    methods
        function self = TestMovieDataSanityCheck(name)
            self = self@TestCase(name);
        end

        %% Set up and tear down methods
        function setUp(self)
            channelPath = @(x) fullfile(self.moviePath, ['Channel' num2str(x)]);
            channels(self.nChan, 1) = Channel();
            for i = 1:self.nChan
                channels(i) = TestHelperMovieObject.setUpChannel(channelPath(i),...
                    self.imSize,self.nFrames);
            end
            self.movie=TestHelperMovieObject.setUpMovie(self.moviePath,channels);
        end
        
        function tearDown(self)
            delete(self.movie);
            rmdir(self.moviePath,'s');
        end
        
        function sanityCheck(self)
            self.movie.sanityCheck()
        end
        
        function newmovie = reload(self)
             newmovie = MovieData.load(self.movie.getFullPath);
        end
        
        %% Regular sanityCheck()            
        function testSanityCheck(self)
            self.movie.sanityCheck();
            assertEqual(self.movie.channels_.owner_,self.movie);
            assertEqual(self.movie.imSize_,self.imSize);
            assertEqual(self.movie.nFrames_,self.nFrames);
            assertEqual(numel(self.movie.channels_),self.nChan);
        end
         
        %% Invalid sizeT tests
        function testInvalidNumberFrames(self)            
            TestHelperMovieObject.setUpChannel(self.movie.channels_(1).channelPath_,...
                    self.imSize(end:-1:1),self.nFrames+1);
            assertExceptionThrown(@() self.sanityCheck, 'MovieData:sanityCheck:nFrames');
        end
        
        function testReloadInvalidNumberFrames(self)
            self.sanityCheck();
            for i = 1: self.nChan
                TestHelperMovieObject.setUpChannel(self.movie.channels_(i).channelPath_,...
                    self.imSize,self.nFrames+1);
            end
            assertExceptionThrown(@() self.reload, 'MovieData:sanityCheck:nFrames');
        end

        %% Invalid sizeX, sizeY tests
        function testInvalidImSize(self)
            TestHelperMovieObject.setUpChannel(self.movie.channels_(1).channelPath_,...
                    self.imSize(end:-1:1),self.nFrames);
            assertExceptionThrown(@() self.sanityCheck, 'MovieData:sanityCheck:imSize');
        end
        
                
        function testReloadInvalidImSize(self)
            self.sanityCheck();
            for i = 1: self.nChan
                TestHelperMovieObject.setUpChannel(self.movie.channels_(i).channelPath_,...
                    self.imSize(end:-1:1),self.nFrames);
            end
            assertExceptionThrown(@() self.reload, 'MovieData:sanityCheck:imSize');
        end
    end
end
