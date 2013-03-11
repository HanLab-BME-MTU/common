% Integration tests for the imSpeckleTracking utility function
%
% Require MATLAB xUnit Test Framework to be installed
% http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework

classdef TestImSpeckleTracking < TestCase
    
    properties
        stack
        imsize = [512 512]
        depth = 2
        x0 = [256 256]
        v0
        minCorL = 20
        maxCorL = 20
        mode = 'fast'
        maxSpd = 10
        verbose = false
        v
    end
    
    
    methods
        function self = TestImSpeckleTracking(name)
            self = self@TestCase(name);
        end
        
        function trackGaussianSpots(self)
            % Generate a imsize x depth stack with 2D Gaussian spots 
            % centered at x0 moving with a velocity v0 (in pixels/frame)
            
            % Initialize Gaussian standard deviatiaion and amplitude
            sigma = 3;
            A = 1e4;
            
            % Generate the stack
            self.stack = zeros([self.imsize self.depth]);
            for i = 1 : self.depth
                self.stack(:,:,i) = simGaussianSpots(...
                    self.imsize(1), self.imsize(2), sigma,...
                    'x' ,self.x0(:,1) + (i - 1) * self.v0(:,1),...
                    'y', self.x0(:,2) + (i - 1) * self.v0(:,2), 'A',A);
            end
            
            % Compute the displacement of templates centered at x0
            self.v = imSpeckleTracking(self.stack, self.x0, self.minCorL,...
                self.maxCorL, 'maxSpd', self.maxSpd, 'mode', self.mode,...
                'verbose', self.verbose);
        end
        
        % Simple displacement tests
        function testZeroDisplacement(self)
            % Test the displacement of a single 2D Gaussian spot
            
            self.v0 = [0 0];
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end
        
        function testUnitXDisplacement(self)
            % Test the displacement of a single 2D Gaussian spot
            
            self.v0 = [1 0];
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end
        
        function testUnitYDisplacement(self)
            % Test the displacement of a single 2D Gaussian spot
            
            self.v0 = [0 1];
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end
        
        function testUnitXYDisplacement(self)
            % Test the displacement of a single 2D Gaussian spot
            
            self.v0 = [1 1];
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end
        
        % Max speed tests
        function testMaxSpeedPassing(self)
            % Test the displacement of a single 2D Gaussian spot
            
            self.v0 = [.95*self.maxSpd .95*self.maxSpd];
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end
        
        function testMaxSpeedFailing(self)
            % Test the displacement of a single 2D Gaussian spot
            
            self.v0 = [1.1*self.maxSpd 1.1*self.maxSpd];
            self.trackGaussianSpots();
            assertTrue(all(isnan(self.v)));
        end
        
        function testFloatXYDisplacement(self)
            % Test the displacement of a single 2D Gaussian spot
            
            self.v0 = [1.5 1.5];
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end
        
        function testMultipleGaussianSpots(self)
            % Test the displacement of a multiple 2D Gaussian spot
            
            nPoints = 4;
            dl = 1/(nPoints+1)*(1:nPoints)';
            self.x0 = horzcat(self.imsize(1) * dl, self.imsize(2) * dl);
            self.v0 = repmat([1 1], nPoints, 1);
            
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end
        
        function testStackDepth(self)
            % Test the depth integration of the cross-correlation score
            
            self.v0 = [1 1];
            self.depth=6;
            self.trackGaussianSpots();
            assertElementsAlmostEqual(self.v, self.v0 ,'relative',1e-3);
        end        
    end
end