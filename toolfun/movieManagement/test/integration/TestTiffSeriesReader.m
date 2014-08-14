classdef TestTiffSeriesReader <  TestCase
    
    properties
        path
        reader
        sizeX = 10
        sizeY = 20
        sizeC = 1
        sizeT = 1
        sizeZ = 1
    end
    
    methods
        function self = TestTiffSeriesReader(name)
            self = self@TestCase(name);
        end
        
        %% Set up and tear down methods
        function setUp(self)
            if isunix,
                self.path = '/tmp/testTiffSeries';
            else
                self.path = 'C:\testTiffSeries';
            end
            mkdir(self.path);
        end
        
        function tearDown(self)
            delete(self.reader);
            if exist(self.path, 'dir') == 7,
                rmdir(self.path, 's');
            end
        end
        
        function checkDimensions(self)
            for c = 1 : self.sizeC
                assertEqual(self.reader.getSizeX(1), self.sizeX);
                assertEqual(self.reader.getSizeY(1), self.sizeY);
                assertEqual(self.reader.getSizeZ(1), self.sizeZ);
                assertEqual(self.reader.getSizeC(1), self.sizeC);
                assertEqual(self.reader.getSizeT(1), self.sizeT);
            end
        end
            
        %% Test data formats
        function testIndividualTiffFiles(self)
            self.sizeT = 5;
            I = ones(self.sizeY, self.sizeX, 'uint8');
            for t = 1 : self.sizeT
                imwrite(t * I, fullfile(self.path, ['test' num2str(t) '.tif']));
            end
            self.reader = TiffSeriesReader({self.path});
            
            self.checkDimensions();
            for t = 1 : self.sizeT
                assertEqual(self.reader.loadImage(1, t, 1), t * I);
            end
        end
        
        function testSingleMultiPageTiff(self)
            self.sizeT = 5;
            I = ones(self.sizeY, self.sizeX, 'uint8');
            imPath = fullfile(self.path, 'test.tif');
            imwrite(I, imPath);
            for i = 2 : self.sizeT
                imwrite(i * I, imPath, 'write', 'append');
            end
            self.reader = TiffSeriesReader({self.path});
            
            self.checkDimensions();
            for t = 1 : self.sizeT
                assertEqual(self.reader.loadImage(1, t, 1), t * I);
            end
        end
        
        function testMultipleMultiPageTiff(self)
            self.sizeT = 5;
            self.sizeZ = 3;
            for t = 1 : self.sizeT
                imPath = fullfile(self.path, ['test' num2str(t) '.tif']);
                I = ones(self.sizeY, self.sizeX, 'uint8');
                imwrite(t * I, imPath);
                for z = 2 : self.sizeZ
                    imwrite(t * I, imPath, 'write', 'append');
                end
            end
            self.reader = TiffSeriesReader({self.path});
            
            self.checkDimensions();
            for t = 1 : self.sizeT
                assertEqual(self.reader.loadImage(1, t, 1), t * I);
            end
        end
        
        %% Test pixel types
        function testMultiChannel(self)
            I = ones(self.sizeY, self.sizeX, 'uint8');
            self.sizeC = 4;
            chPath = cell(self.sizeC, 1);
            for c = 1 : self.sizeC
                chPath{c} = fullfile(self.path, ['ch' num2str(c)]);
                mkdir(chPath{c});
                imwrite(c * I, fullfile(chPath{c}, 'test.tif'));
            end
            self.reader = TiffSeriesReader(chPath);
            
            self.checkDimensions();
            for c = 1 : self.sizeC
                assertEqual(self.reader.loadImage(c, 1, 1), c * I);
            end
        end
        %% Test pixel types
        function testUINT8(self)
            I = ones(self.sizeY, self.sizeX, 'uint8');
            imwrite(I, fullfile(self.path, 'test.tif'));
            self.reader = TiffSeriesReader({self.path});
            
            assertEqual(self.reader.getBitDepth(1), 8);
            assertEqual(self.reader.loadImage(1, 1, 1), I);
        end
        
        function testUINT16(self)
            I = ones(self.sizeY, self.sizeX, 'uint16');
            imwrite(I, fullfile(self.path, 'test.tif'));
            self.reader = TiffSeriesReader({self.path});
            
            assertEqual(self.reader.getBitDepth(1), 16);
            assertEqual(self.reader.loadImage(1, 1, 1), I);
        end
    end
end
