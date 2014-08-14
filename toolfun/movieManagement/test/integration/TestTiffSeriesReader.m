classdef TestTiffSeriesReader <  TestCase
    
    properties
        path
        reader
        sizeX = 10
        sizeY = 20
        sizeC = 1
        sizeT = 1
        sizeZ = 1
        imClass = 'uint8'
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
        
        function checkLoadImage(self)
            for c = 1 : self.sizeC
                for t = 1 : self.sizeT
                    for z = 1 : self.sizeZ
                        assertEqual(self.reader.loadImage(c, t, z),...
                            self.getPlane(c, t, z));
                    end
                end
            end
        end
        
        function I = getPlane(self, c, t, z)
            index = sub2ind([self.sizeC self.sizeT self.sizeZ], c, t, z);
            I = index * ones(self.sizeY, self.sizeX, self.imClass);
        end
        
        %% Test data formats
        function testIndividualTiffFiles(self)
            self.sizeT = 5;
            for t = 1 : self.sizeT
                imwrite(self.getPlane(1, t, 1),...
                    fullfile(self.path, ['test' num2str(t) '.tif']));
            end
            self.reader = TiffSeriesReader({self.path});
            
            self.checkDimensions();
            self.checkLoadImage();
            assertFalse(self.reader.isSingleMultiPageTiff(1));
        end
        
        function testSingleMultiPageTiff(self)
            self.sizeT = 5;
            imPath = fullfile(self.path, 'test.tif');
            imwrite(self.getPlane(1, 1, 1), imPath);
            for t = 2 : self.sizeT
                imwrite(self.getPlane(1, t, 1), imPath, 'write', 'append');
            end
            self.reader = TiffSeriesReader({self.path});
            
            self.checkDimensions();
            self.checkLoadImage();
            assertTrue(self.reader.isSingleMultiPageTiff(1));
        end
        
        function testMultipleMultiPageTiff(self)
            self.sizeT = 5;
            self.sizeZ = 3;
            for t = 1 : self.sizeT
                imPath = fullfile(self.path, ['test' num2str(t) '.tif']);
                imwrite(self.getPlane(1, t, 1), imPath);
                for z = 2 : self.sizeZ
                    imwrite(self.getPlane(1, t, z),...
                        imPath, 'write', 'append');
                end
            end
            self.reader = TiffSeriesReader({self.path});
            
            self.checkDimensions();
            self.checkLoadImage();
            assertFalse(self.reader.isSingleMultiPageTiff(1));
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
