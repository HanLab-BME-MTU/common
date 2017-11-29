classdef TestBootStrapMean < TestCase
    %TestBootStrapMean Tests the function bootStrapMean
    
    properties
        % Number of which to find mean
        numbers
        % Alpha value
        alpha
        % Number of bootstraps to do
        nBoot
    end
    
    methods
        function self = TestBootStrapMean(name)
            self@TestCase(name);
        end
        function setUp(self)
            self.numbers = rand(100,100);
            self.alpha = 0.05;
            self.nBoot = 1000;
        end
        function tearDown(self)
            % Delete parallel pool on tear down
            delete(gcp('nocreate'))
        end
        function TestWithoutParallel(self)
            % Make sure no parallel pool exists
            delete(gcp('nocreate'));
            [conf,meanS] = bootStrapMean(self.numbers, self.alpha, self.nBoot);
            meanO = mean(self.numbers);
            assert(all(meanS > conf(1,:)));
            assert(all(meanS < conf(2,:)));
            assert(all(meanO > conf(1,:)));
            assert(all(meanO < conf(2,:)));
        end
        function TestWithParallel(self)
            % Make sure parallel pool exists
            if(isempty(gcp('nocreate')))
                parpool(3);
            end
            [conf,meanS] = bootStrapMean(self.numbers, self.alpha, self.nBoot);
            meanO = mean(self.numbers);
            assert(all(meanS > conf(1,:)));
            assert(all(meanS < conf(2,:)));
            assert(all(meanO > conf(1,:)));
            assert(all(meanO < conf(2,:)));
        end
    end
    
end

