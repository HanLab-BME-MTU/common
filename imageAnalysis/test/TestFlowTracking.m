% Integration tests for the imSpeckleTracking utility function
%
% Require MATLAB xUnit Test Framework to be installed
% http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework

classdef TestFlowTracking < TestCase
    
    properties
        imsize = [512 512];
        stack
        f = @imSpeckleTracking
    end
    
    
    methods
        function self = TestFlowTracking(name)
            self = self@TestCase(name);
        end
        
        function testSingleGaussianSpot(self)
        % Test the displacement of a single 2D Gaussian spot
        
            v0 = [0 0; 1 0; 0 1; 1 1; .5 .5; 5.5 5.5; 9.5 9.5];
            for i = 1 : size(v0,1)
                self.setUpGaussianSpot(v0(i,:));
                v = self.f(self.stack, self.imsize/2, 20, 20);
                assertElementsAlmostEqual(v, v0(i,:) ,'relative',1e-3);
            end
        end
        
        function testMaxFlowSpeed(self)
            % Test the maximum flow speed parameter

            % Set up a single 2D Gaussian spot moving by [v0 v0]
            v0 = 5;
            vth = repmat(v0, 1, 2);
            self.setUpGaussianSpot(vth);

            % Test maximum speeds larger than v0
            for maxSpd = v0+1 : 2 : 20;
                v = self.f(self.stack, self.imsize/2, 20, 20,...
                    'maxSpd', maxSpd);
                assertElementsAlmostEqual(v, vth, 'relative', 1e-4);
            end
            
            % Test maximum speeds smaller than v0
            for maxSpd = 1 : .5 : v0 - .5
                v = self.f(self.stack, self.imsize/2, 20, 20,...
                    'maxSpd', maxSpd);
                assertTrue(all(isnan(v)));
            end
        end

    
        function setUpGaussianSpot(self, v) 
            % Generate a n x m x2 stack with a single 2D Gaussian spot
            % moving with a velocity v
            
            assert(numel(v) == 2, 'Velocity must be a 2D vector');
            
            sigma = 3;
            A = 1e4;
            self.stack = zeros([self.imsize 2]);
            self.stack(:,:,1) = simGaussianSpots(self.imsize(1), self.imsize(2), sigma,...
                'x' ,self.imsize(1)/2, 'y', self.imsize(2)/2, 'A',A);
            self.stack(:,:,2) = simGaussianSpots(self.imsize(1), self.imsize(2), sigma,...
                'x' ,self.imsize(1)/2+v(1), 'y', self.imsize(2)/2+v(2), 'A',A);
        end
    end

end