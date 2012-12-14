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
                assertElementsAlmostEqual(v, v0(i,:) ,'relative',1e-3,...
                    ['Displacement: ' num2str(v0(i,:))]);
            end
        end
        
        
        function testStackDepth(self)
            % Test the depth integration of the cross-correlation score
            
            v0 = [1 1];
            for depth = 2 : 15
                self.setUpGaussianSpot(v0, depth);
                v = self.f(self.stack, self.imsize/2, 20, 20);
                assertElementsAlmostEqual(v, v0 ,'relative',1e-2,...
                    ['Stack depth: ' num2str(depth)]);
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
                assertElementsAlmostEqual(v, vth, 'relative', 1e-4,...
                    ['Maximum speed: ' num2str(maxSpd)]);
            end
            
            % Test maximum speeds smaller than v0
            for maxSpd = 1 : .5 : v0 - .5
                v = self.f(self.stack, self.imsize/2, 20, 20,...
                    'maxSpd', maxSpd);
                assertTrue(all(isnan(v)),...
                    ['Maximum speed: ' num2str(maxSpd)]);
            end
        end
        
        
        function setUpGaussianSpot(self, v, varargin)
            % Generate a n x m x2 stack with a single 2D Gaussian spot
            % moving with a velocity v (in pixels/frame)
            
            % Check input
            ip = inputParser;
            ip.addRequired('v', @(x) isvector(x) && numel(v) == 2);
            ip.addOptional('depth', 2, @isscalar);
            ip.parse(v, varargin{:});
            depth = ip.Results.depth;
            
            % Initialize Gaussian standard deviatiaion and amplitude
            sigma = 3;
            A = 1e4;
            
            % Generate the stack
            self.stack = zeros([self.imsize depth]);
            for i = 1 : depth
                self.stack(:,:,i) = simGaussianSpots(...
                    self.imsize(1), self.imsize(2), sigma,...
                    'x' ,self.imsize(1)/2 + (i - 1) * v(1),...
                    'y', self.imsize(2)/2 + (i - 1) * v(2), 'A',A);
            end
        end
    end
end