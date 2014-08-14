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
    end
end
