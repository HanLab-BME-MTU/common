classdef TestTiffSeriesReader <  TestCase
    
    properties
        path
        reader
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
        
        %% Test data formats
        function testIndividualTiffFiles(self)
            I = ones(10, 10, 'uint8');
            for t = 1 :5
                imwrite(t * I, fullfile(self.path, ['test' num2str(t) '.tif']));
            end
            self.reader = TiffSeriesReader({self.path});
            
            assertEqual(self.reader.getSizeX(1), 10);
            assertEqual(self.reader.getSizeY(1), 10);
            assertEqual(self.reader.getSizeZ(1), 1);
            assertEqual(self.reader.getSizeC(1), 1);
            assertEqual(self.reader.getSizeT(1), 5);
            
            for t = 1 : 5
                assertEqual(self.reader.loadImage(1, t, 1), t * I);
            end
        end
        
        function testSingleMultiPageTiff(self)
            I = ones(10, 10, 'uint8');
            imPath = fullfile(self.path, 'test.tif');
            imwrite(I, imPath);
            for i = 2 : 5
                imwrite(i * I, imPath, 'write', 'append');
            end
            self.reader = TiffSeriesReader({self.path});
            
            assertEqual(self.reader.getSizeX(1), 10);
            assertEqual(self.reader.getSizeY(1), 10);
            assertEqual(self.reader.getSizeZ(1), 1);
            assertEqual(self.reader.getSizeC(1), 1);
            assertEqual(self.reader.getSizeT(1), 5);
            
            for t = 1 : 5
                assertEqual(self.reader.loadImage(1, t, 1), t * I);
            end
        end
        
        function testMultipleMultiPageTiff(self)
            for t = 1 : 5
                imPath = fullfile(self.path, ['test' num2str(t) '.tif']);
                I = ones(10, 10, 'uint8');
                imwrite(t * I, imPath);
                imwrite(t * I, imPath, 'write', 'append');
            end
            self.reader = TiffSeriesReader({self.path});
            
            assertEqual(self.reader.getSizeX(1), 10);
            assertEqual(self.reader.getSizeY(1), 10);
            assertEqual(self.reader.getSizeZ(1), 2);
            assertEqual(self.reader.getSizeC(1), 1);
            assertEqual(self.reader.getSizeT(1), 5);
            
            for t = 1 : 5
                assertEqual(self.reader.loadImage(1, t, 1), t * I);
            end
        end
        
        %% Test pixel types
        function testMultiChannel(self)
            I = ones(10, 10, 'uint8');
            chPath = cell(2, 1);
            for c = 1 : 2
                chPath{c} = fullfile(self.path, ['ch' num2str(c)]);
                mkdir(chPath{c});
                imwrite(c * I, fullfile(chPath{c}, 'test.tif'));
            end
            self.reader = TiffSeriesReader(chPath);
            
            assertEqual(self.reader.getSizeX(1), 10);
            assertEqual(self.reader.getSizeY(1), 10);
            assertEqual(self.reader.getSizeZ(1), 1);
            assertEqual(self.reader.getSizeC(1), 2);
            assertEqual(self.reader.getSizeT(1), 1);
            
            for c = 1 : 2
                assertEqual(self.reader.loadImage(c, 1, 1), c * I);
            end
        end
        %% Test pixel types
        function testUINT8(self)
            I = ones(10, 10, 'uint8');
            imwrite(I, fullfile(self.path, 'test.tif'));
            self.reader = TiffSeriesReader({self.path});
            
            assertEqual(self.reader.getBitDepth(1), 8);
            assertEqual(self.reader.loadImage(1, 1, 1), I);
        end
        
        function testUINT16(self)
            I = ones(10, 10, 'uint16');
            imwrite(I, fullfile(self.path, 'test.tif'));
            self.reader = TiffSeriesReader({self.path});
            
            assertEqual(self.reader.getBitDepth(1), 16);
            assertEqual(self.reader.loadImage(1, 1, 1), I);
        end
    end
end
