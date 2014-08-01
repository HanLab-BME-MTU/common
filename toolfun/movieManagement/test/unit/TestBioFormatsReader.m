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
        
        function testConstructorSeries(self)
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
    end
end
