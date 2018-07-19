classdef TestInterpftExtrema < TestCase
    %TESTINTERPFTEXTREMA Test function interpft_extrema
    
    properties
        seed
        figures
    end
    
    methods
        function self = TestInterpftExtrema(varargin)
            self = self@TestCase(varargin{:});
        end
        function setUp(self)
            self.seed = 7.365476979391746e+05;
            rng(self.seed);
            self.figures = [];
        end
        function tearDown(self)
            close(self.figures);
        end
        function testOneDimensionalExample(self)
          r = rand(7,1);
          h = figure;
          plot((0:length(r)-1)/length(r)*2*pi,r,'ko');
          hold on;
          plot((0:199)/200*2*pi,interpft(r,200),'k');
          interpft_extrema(r);
          hold off;
          self.figures(end+1) = h;
        end
        function testTwoDimensionalExample(self)
          r = rand(11,3);
          h = figure;
          plot((0:size(r,1)-1)/size(r,1)*2*pi,r,'ko');
          hold on;
          plot((0:199)/200*2*pi,interpft(r,200),'k');
          interpft_extrema(r);
          hold off;
          self.figures(end+1) = h;
        end
        function testThreeDimensionalExample(self)
            r = rand(11,5,4);
            [maxima, minima, maxima_value, minima_value, other, other_value] = interpft_extrema(r);
        end
        function testEvenInput(self)
            r = rand(14,1);
            [maxima, minima, maxima_value, minima_value, other, other_value] = interpft_extrema(r);
        end
        function testDimenions(self)
            r = rand(12,4);
            maxima_first = interpft_extrema(r);
            maxima_second = interpft_extrema(permute(r,[3 2 1]),3);
            assertEqual(permute(maxima_first,[3 2 1]),maxima_second);
        end
        function testSort(self)
            r = rand(17,2);
            [maxima, minima, maxima_value, minima_value, other, other_value] = interpft_extrema(r,1,false);
            [maxima_s, minima_s, maxima_value_s, minima_value_s, other_s, other_value_s] = interpft_extrema(r,1,true);
            assertEqual(nanmax(maxima_value),maxima_value_s(1,:));
        end
        function testSinInput(self)
            for n=5:10
                x = 0:n-1;
                x = x'/n*2*pi;
                [maxima, minima, maxima_value, minima_value, other, other_value] = interpft_extrema(sin(x),1,true);
                assert(abs(maxima_value-1) < 1e-15);
                assert(abs(minima_value+1) < 1e-15);
                assertEqual(sin(maxima),1);
                assertEqual(sin(minima),-1);
            end
        end
        function testCosInput(self)
            for n=5:10
                x = 0:n-1;
                x = x'/n*2*pi;
                % Relax the tolerance since the maximum occurs at the
                % boundary
                [maxima, minima, maxima_value, minima_value, other, other_value] = interpft_extrema(cos(x),1,true,1e-10);
                assert(abs(maxima_value-1) < 1e-15);
                assert(abs(minima_value+1) < 1e-15);
                assertEqual(cos(maxima),1);
                assertEqual(cos(minima),-1);
            end
        end
        function testConstant(self)
            maxima = interpft_extrema(ones(5,3),1,true);
            assert(isempty(maxima));
        end
        function testMixed(self)
            r = rand(16,6);
            r(:,4) = 5;
            maxima = interpft_extrema(r);
            assert(all(isnan(maxima(:,4))));
        end
    end

    
end

