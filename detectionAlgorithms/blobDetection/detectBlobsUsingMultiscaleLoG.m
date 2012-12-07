function [ imBlobLocations, varargout ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
% Detects blobs as local maxima of the Multiscale LoG filtered image
% 
%   [ imBlobLocations ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
%   [ imBlobLocations, imMultiscaleLoGResponse ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
%   [ imBlobLocations, imMultiscaleLoGResponse, imPixelScaleMap ] = detectBlobsUsingMultiscaleLoG( im, blobDiameterRange, varargin )
% 
% The input intensity image is first filtered with a Laplacian of Gaussian
% (LoG) Filter accross multiple scales/sigmas.
% The seed points are then detected as local maxima in the scale space.
%
% References:
%   
% Byun, J., M. R. Verardo, et al. (2006). 
% "Automated tool for the detection of cell nuclei in digital microscopic 
% images: application to retinal images." Mol Vis 12: 949-960.
%
% Author: Deepak Roy Chittajallu
%

    p = inputParser;    
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.addRequired( 'blobDiameterRange', @(x) (numel(x) == 2) );    
    p.parse( im, blobDiameterRange );    
    
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'minBlobDistance', [], @(x) isscalar(x) );
    p.addParamValue( 'numLoGScales', 15, @(x) isscalar(x) );
    p.addParamValue( 'flagBrightBlobs', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, blobDiameterRange, varargin{:} ); 
    
    minBlobDistance = p.Results.minBlobDistance;
    spacing = p.Results.spacing;    
    numLoGScales = p.Results.numLoGScales;
    flagBrightBlobs = p.Results.flagBrightBlobs;
    flagDebugMode = p.Results.debugMode;

    % Compute LoG at a series of sigmas
    sigmaLogRange = log2( (0.5 * sort(blobDiameterRange) / sqrt(ndims(im))) );
    sigmaValues = 2.^linspace( sigmaLogRange(1), sigmaLogRange(2), numLoGScales);

    [ imMultiscaleLoGResponse, pixelScaleMap ] = filterMultiscaleLoGND( im, sigmaValues, ...
                                                                        'spacing', spacing, ...
                                                                        'debugMode', flagDebugMode );

    if flagBrightBlobs
        imMultiscaleLoGResponse = -1 * imMultiscaleLoGResponse;
    end
    
    % locate local intensity maxima in gaussian blurred image
    if isempty(minBlobDistance)
        MaximaSuppressionSize = round(  min( blobDiameterRange ) ./ spacing );
    else
        MaximaSuppressionSize = round( minBlobDistance ./ spacing );    
    end
    
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    MaximaSuppressionSize( MaximaSuppressionSize < 3 ) = 3; 
    
    switch ndims( im ) 
        
        case 2 
            
            imLocalMax = locmax2d(imMultiscaleLoGResponse, MaximaSuppressionSize, 1);            
            
            if flagDebugMode
                
                seedInd = find( imLocalMax > 0 );
                [cy, cx] = ind2sub( size(imLocalMax), seedInd );                
                theta = 0:0.1:(2*pi+0.1);
                
                cx = cx(:,ones(size(theta)));
                cy = cy(:,ones(size(theta)));
                rad = (sigmaValues(pixelScaleMap(seedInd)))' * sqrt(ndims(im));
                rad = rad(:,ones(size(theta)));
                
                theta = theta(ones(size(cx,1),1),:);                
                X = cx + cos(theta).* rad;
                Y = cy + sin(theta).* rad;

                imseriesmaskshow( imMultiscaleLoGResponse, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Multiscale LoG Filter' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
                imseriesmaskshow( im, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Seed Points Overlayed on Input Image' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
                figure, histogram( pixelScaleMap(seedInd) );
                title( 'Histogram of the cell scales (diameter) found in the image' );
                
            end
            
        case 3
            
            imLocalMax = locmax3d(imMultiscaleLoGResponse, MaximaSuppressionSize);           
            
            if flagDebugMode
                imseriesmaskshow( imMultiscaleLoGResponse, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Multiscale LoG Filter' );
                
                imseriesmaskshow( im, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Seed Points Overlayed on Input Image' );
            end
            
    end

    % detect local intensity maxima as cell seed points
    imBlobLocations = imLocalMax;   
    if nargout > 1
        varargout{1} = imMultiscaleLoGResponse;
    end 
    
    if nargout > 2
        % return estimated scale or radius for each detected blob
        varargout{2} = sigmaValues( pixelScaleMap ) * sqrt(ndims(im)); 
    end
    
end