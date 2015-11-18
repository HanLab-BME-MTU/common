classdef TestRobustMean < TestCase
    %TESTROBUSTMEAN Test the robustMean functions
    
    properties
        distribution
    end
    
    methods
        function self = TestRobustMean(name)
            self = self@TestCase(name);
        end
        function setUp(self)
            % Seed for consistency
            rng(39595353);
            self.distribution = randn(100);
        end
        function testVector(self)
            v = self.distribution(:)';
            [initRobustMean,initRobustStd] = robustMean(v);
            initMean = mean(v);
            initStd = std(v);
            assert(abs(initRobustMean -initMean) < 1e-1);
            assert(abs(initRobustStd - initStd) < 1e-1);
            for i=1:10;
                nv = [v ones(1,i)*100];
                m = mean(nv);
                s = std(nv);
                [rm,rs] = robustMean(nv);
                assert(abs(rm-initRobustMean) < abs(m-initMean));
                assert(abs(rs-initRobustStd) < abs(s-initStd));
            end
        end
        function testMatrix(self)
            v = self.distribution;
            [initRobustMean,initRobustStd] = robustMean(v);
            initMean = mean(v);
            initStd = std(v);
            assert(abs(mean(initRobustMean -initMean)) < 1e-1);
            assert(abs(mean(initRobustStd - initStd)) < 1e-1);
            for i=1:10;
                nv = [v ; ones(i,size(v,2))*100];
                m = mean(nv);
                s = std(nv);
                [rm,rs] = robustMean(nv);
                assert(all(abs(rm-initRobustMean) < abs(m-initMean)));
                assert(all(abs(rs-initRobustStd) < abs(s-initStd)));
            end
        end
    end
    
end

