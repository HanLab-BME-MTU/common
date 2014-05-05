classdef TestBFMovieData < TestMovieData & TestCase
    
    properties
        fakename = 'test.fake';
        lociToolsPath
    end
    
    methods
        function self = TestBFMovieData(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            self.setUp@TestMovieData();
            
            % Get path to Bio-Formats JAR (assuming it is in Matlab path)
            self.lociToolsPath = which('bioformats_package.jar');
            assert(~isempty(self.lociToolsPath));
            
            % Remove Bio-Formats JAR from dynamic class path
            if ismember(self.lociToolsPath,javaclasspath('-dynamic'))
                javarmpath(self.lociToolsPath);
            end
            
            bfCheckJavaPath;
            r = loci.formats.in.FakeReader();
            self.imSize = [r.DEFAULT_SIZE_Y r.DEFAULT_SIZE_X];
            self.nChan = r.DEFAULT_SIZE_C;
            self.nFrames = r.DEFAULT_SIZE_T;
        end
        
        function tearDown(self)
            self.tearDown@TestMovieData();
            
            % Remove Bio-Formats JAR from dynamic class path
            if ismember(self.lociToolsPath,javaclasspath('-dynamic'))
                javarmpath(self.lociToolsPath);
            end
        end
        
        function filename = createFakeFile(self)
            filename = fullfile(self.path, self.fakename);
            fid = fopen(filename, 'w');
            fclose(fid);
        end
        
        function filename = createFakeFileCompanion(self, content)
            filename = fullfile(self.path, [self.fakename '.ini']);
            fid = fopen(filename, 'w');
            fwrite(fid, content);
            fclose(fid);
        end
        
        function setUpMovie(self)
            filename = self.createFakeFile();
            self.movie = MovieData(filename);
        end
        
        function checkChannelPaths(self)
            for i = 1 : self.nChan
                assertEqual(self.movie.getChannel(i).channelPath_,...
                    fullfile(self.path, self.fakename))
            end
        end
        
        %% Constructor
        
        function testConstructor(self)
            filename = self.createFakeFile();
            self.movie = MovieData(filename);
            self.checkChannelPaths();
        end
        
        function testConstructorMetadata(self)
            filename = self.createFakeFile();
            self.createFakeFileCompanion('physicalSizeX=1');
            self.movie = MovieData(filename);
            self.checkChannelPaths();
        end
        
        function testConstructorMultiSeries(self)
            self.fakename = 'test&series=2.fake';
            filename = self.createFakeFile();
            movies = MovieData(filename);
            assertEqual(numel(movies), 2);
            assertFalse(isequal(movies(1).getReader().formatReader,...
                movies(2).getReader().formatReader));
        end
        
        function testConstructorMultiSeriesReuseReader(self)
            self.fakename = 'test&series=2.fake';
            filename = self.createFakeFile();
            movies = MovieData(filename, 'reuseReader', true);
            assertEqual(numel(movies), 2);
            assertEqual(movies(1).getReader().formatReader,...
                movies(2).getReader().formatReader);
        end
        
        %% Typecasting tests
        function checkPixelType(self, classname)
            if strcmp(classname, 'single'),
                pixelsType = 'float';
            else
                pixelsType = classname;
            end
            self.fakename = ['test&pixelType=' pixelsType '.fake'];
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
        function testSizeX(self)
            self.fakename = 'test&sizeX=100.fake';
            self.imSize(2) = 100;
            self.setUpMovie()
            self.checkDimensions();
        end
        
        function testSizeY(self)
            self.fakename = 'test&sizeY=100.fake';
            self.imSize(1) = 100;
            self.setUpMovie()
            self.checkDimensions();
        end
        
        function testSizeZ(self)
            self.fakename = 'test&sizeZ=256.fake';
            self.zSize = 256;
            self.setUpMovie()
            self.checkDimensions();
        end
        
        function testSizeC(self)
            self.fakename = 'test&sizeC=4.fake';
            self.nChan = 4;
            self.setUpMovie()
            self.checkDimensions();
        end
        
        function testSizeT(self)
            self.fakename = 'test&sizeT=256.fake';
            self.nFrames = 256;
            self.setUpMovie()
            self.checkDimensions();
        end
        
        %% ROI tests
        function testAddROIMultiSeries(self)
            nMovies = 3;
            self.fakename = ['test&series=' num2str(nMovies) '.fake'];
            self.setUpMovie();
            assertEqual(numel(self.movie), nMovies);
            for i = 1 : nMovies
                self.movie(i).addROI('','');
                roi = self.movie(i).getROI(1);
                assertEqual(roi.getChannel(1), self.movie(i).getChannel(1));
                assertEqual(roi.getSeries(), self.movie(i).getSeries());
            end
        end
    end
end
