function [res, theta, nms, pixelScaleMap] = multiscaleSteerableDetector3D( im, filterType, sigmaValues, varargin )
% Performs multiscale steerable filtering in 3D 
%
% A steerable filter is applied in multiple scales and results are combined
% by selecting the per-pixel maximum response across scales.
% The resulting scale map is used to produce an orientation map.
%
% Inputs:
% 
%                    im : input 3D volume
%             
%            filterType : type of steerable filter
% 
%                         1: curve detector
%                         2: surface detector
%                     
%           sigmaValues : An array of scales (standard deviation of the 
%                         gaussian kernel) on which the steerable filter 
%                         will be applied.
%                      
%     zAnisotropyFactor : zSpacing / xySpacing (optional argument)
%                         Default: Assumes isotropic volume  
%                      
% Outputs:
% 
%              res : response of the multi-scale steerable filter
%            theta : orientation map
%              nms : non-maximal suppression response
%    pixelScaleMap : pixel map of scales for which response is maximal
% 
% See also steerableDetector3D.m, testSteerableDetector3D.m
% 
% Author: Deepak Roy Chittajallu with help from Francois Aguet 
% This function is based on steerableDetector3D.m written by Francois Aguet
%

    p = inputParser;
    p.addRequired( 'im', @(x) ( isnumeric(x) && ndims(x) == 3 ) );
    p.addRequired( 'filterType', @(x) (ismember(x, [1,2])) );
    p.addRequired( 'sigmaValues', @(x) isnumeric(x) );
    p.addOptional( 'zAnisotropyFactor', 1.0, @(x) isscalar(x) );
    p.parse( im, filterType, sigmaValues, varargin{:} );
    
    imsize = size(im);
    zAnisotropyFactor = p.Results.zAnisotropyFactor;

    fprintf( '\nRunning steerable detector at multiple scales on %d x %d x %d sized volume ...\n', imsize(2), imsize(1), imsize(3) );
    for i = 1:numel( sigmaValues )

        fprintf( '\n\t%d/%d: Trying sigma value of %.2f ... ', i, numel( sigmaValues ), sigmaValues(i) );   

        tic
        [curRes, cutTheta, curNms] = steerableDetector3D(im, filterType, sigmaValues(i), zAnisotropyFactor);
        timeElapsed = toc;

        fprintf( 'It took %.2f seconds\n', timeElapsed );   
                
        %curRes = sigmaValues(i)^2 * curRes; % scale normalization
        
        if i == 1        
            res = curRes;
            nms = curNms;
            theta = cutTheta;
            pixelScaleMap = ones( size(res) );
        else
            imBetterMask = curRes > res;
            res(imBetterMask) = curRes(imBetterMask);
            nms(imBetterMask) = curNms(imBetterMask);
            theta.x1(imBetterMask) = cutTheta.x1(imBetterMask);
            theta.x2(imBetterMask) = cutTheta.x2(imBetterMask);
            theta.x3(imBetterMask) = cutTheta.x3(imBetterMask);
            pixelScaleMap(imBetterMask) = i;
        end
    end

end