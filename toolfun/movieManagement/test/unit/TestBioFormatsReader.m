classdef TestBioFormatsReader < TestCase
    
    properties
        id = 'test.fake';
        reader
    end
    
    methods
        function self = TestBioFormatsReader(name)
            self = self@TestCase(name);
        end
        
        %% Test constructor
        function testDefaultConstructor(self)
            self.reader = BioFormatsReader(self.id);
            assertEqual(self.reader.id, self.id);
            assertEqual(self.reader.series, 0);
            assertTrue(isa(self.reader.formatReader,...
                'loci.formats.IFormatReader'));
        end
        
        function testConstructorMultiSeries(self)
            self.id = 'test&series=2.fake';
            self.reader = BioFormatsReader(self.id, 1);
            assertEqual(self.reader.id, self.id);
            assertEqual(self.reader.series, 1);
            assertTrue(isa(self.reader.formatReader,...
                'loci.formats.IFormatReader'));
        end
        
        function testConstructorReader(self)
            r = bfGetReader(self.id);
            self.reader = BioFormatsReader(self.id, 'reader', r);
            assertEqual(self.reader.id, self.id);
            assertEqual(self.reader.series, 0);
            assertEqual(self.reader.formatReader, r);
        end
        
        function testConstructorMultiSeriesReader(self)
            self.id = 'test&series=2.fake';
            r = bfGetReader(self.id);
            self.reader = BioFormatsReader(self.id, 0, 'reader', r);
            self.reader(2) = BioFormatsReader(self.id, 1, 'reader', r);
            assertEqual(self.reader(1).id, self.id);
            assertEqual(self.reader(2).id, self.id);
            assertEqual(self.reader(1).series, 0);
            assertEqual(self.reader(2).series, 1);
            assertEqual(self.reader(1).formatReader, r);
            assertEqual(self.reader(2).formatReader, r);
        end
        
        %% Test getReader
        function testGetReader(self)
            self.reader = BioFormatsReader(self.id);
            assertTrue(isa(self.reader.getReader(),...
                'loci.formats.IFormatReader'));
            assertEqual(char(self.reader.formatReader.getCurrentFile()), self.id);
        end
        
        function testGetReaderInit(self)
            r = bfGetReader(self.id);
            self.reader = BioFormatsReader(self.id, 'reader', r);
            assertEqual(self.reader.getReader(), r);
            assertEqual(char(self.reader.formatReader.getCurrentFile()), self.id);
        end
        
        function testGetReaderMultiSeries(self)
            self.id = 'test&series=2.fake';
            self.reader = BioFormatsReader(self.id, 0);
            self.reader(2) = BioFormatsReader(self.id, 1);
            assertEqual(self.reader(1).formatReader.getSeries(), 0);
            assertEqual(self.reader(2).formatReader.getSeries(), 0);
            assertEqual(self.reader(1).getReader().getSeries(), 0);
            assertEqual(self.reader(2).getReader().getSeries(), 1);
        end
        
        function testGetReaderMultiSeriesReader(self)
            self.id = 'test&series=2.fake';
            r = bfGetReader(self.id);
            self.reader = BioFormatsReader(self.id, 0, 'reader', r);
            self.reader(2) = BioFormatsReader(self.id, 1, 'reader', r);
            assertEqual(char(self.reader(1).formatReader.getCurrentFile()), self.id);
            assertEqual(char(self.reader(2).formatReader.getCurrentFile()), self.id);
            assertEqual(self.reader(1).formatReader.getSeries(), 0);
            assertEqual(self.reader(2).formatReader.getSeries(), 0);
            assertEqual(self.reader(1).getReader().getSeries(), 0);
            assertEqual(self.reader(2).getReader().getSeries(), 1);
        end
        
        function testGetReaderClose(self)
            self.reader = BioFormatsReader(self.id);
            assertTrue(isa(self.reader.getReader(),...
                'loci.formats.IFormatReader'));
            assertEqual(char(self.reader.formatReader.getCurrentFile()), self.id);
            self.reader.formatReader.close();
            assertTrue(isempty(self.reader.formatReader.getCurrentFile()));
            self.reader.getReader();
            assertEqual(char(self.reader.formatReader.getCurrentFile()), self.id);
        end
        
        function testGetReaderMultiSeriesClose(self)
            self.id = 'test&series=2.fake';
            r = bfGetReader(self.id);
            self.reader = BioFormatsReader(self.id, 0, 'reader', r);
            self.reader(2) = BioFormatsReader(self.id, 1, 'reader', r);
            self.reader(2).formatReader.close();
            assertTrue(isempty(self.reader(1).formatReader.getCurrentFile()));
            assertTrue(isempty(self.reader(2).formatReader.getCurrentFile()));
            self.reader(1).getReader();
            assertEqual(char(self.reader(1).formatReader.getCurrentFile()), self.id);
            assertEqual(char(self.reader(2).formatReader.getCurrentFile()), self.id);
        end
    end
end
