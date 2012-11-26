function [ imBlobLocations, varargout ] = detectBlobsUsingMultiscaleLoBG( im, blobDiameterRange, varargin )
% Detects blobs as local maxima of the Multiscale LoBG filtered image
% 
%   [ imBlobLocations ] = detectBlobsUsingMultiscaleLoBG( im, blobDiameterRange, varargin )
%   [ imBlobLocations, imMultiscaleLoBGResponse ] = detectBlobsUsingMultiscaleLoBG( im, blobDiameterRange, varargin )
% 
%   The input intensity image is first filtered with a Laplacian of Bi-Gaussian
%   (LoBG) Filter accross multiple scales/sigmas. The blobs are then detected 
%   as local maxima in the scale space.
%
%   References:
% 
%   Xiao, C., M. Staring, et al. (2012). 
%   "A multiscale bi-Gaussian filter for adjacent curvilinear structures 
%   detection with application to vasculature images." 
%   IEEE Transactions on Image Processing, PP(99): 1-1.
% 
%   See: filterMultiscaleLoBGND.m, filterLoBGND.m
% 
%   Author: Deepak Roy Chittajallu
%

    p = inputParser;    
    p.addRequired( 'im', @(x) (isnumeric(x) && ismember( ndims(x), [2,3] )) );
    p.addRequired( 'blobDiameterRange', @(x) (numel(x) == 2) );    
    p.parse( im, blobDiameterRange );    
    
    p.addParamValue( 'rho', 0.2, @(x) (isscalar(x)) ); % for now read the reference paper to understand what this means
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'numLoGScales', 15, @(x) isscalar(x) );
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, blobDiameterRange, varargin{:} ); 

    spacing = p.Results.spacing;
    numLoGScales = p.Results.numLoGScales;
    flagDebugMode = p.Results.debugMode;
    rho = p.Results.rho;
    
    % compute LoG at a series of sigmas
    sigmaLogRange = log2( 0.5 * sort(blobDiameterRange) );
    sigmaValues = 2.^linspace( sigmaLogRange(1), sigmaLogRange(2), numLoGScales);

    [ imMultiscaleLoBGResponse, pixelScaleMap ] = filterMultiscaleLoBGND( im, sigmaValues, rho, ...
                                                                          'spacing', spacing, ...
                                                                          'debugMode', flagDebugMode );

    imMultiscaleLoBGResponse = -1 * imMultiscaleLoBGResponse;
    
    % locate local intensity maxima in gaussian blurred image
    minCellDiameterImsp = min( blobDiameterRange ) ./ spacing;
    MaximaSuppressionSize = round( minCellDiameterImsp );
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    
    switch ndims( im ) 
        
        case 2 
            
            imLocalMax = locmax2d(imMultiscaleLoBGResponse, MaximaSuppressionSize, 1);            
            imLocalMax = double(imLocalMax > 0);
            
            if flagDebugMode
                
                seedInd = find( imLocalMax > 0 );
                [cy, cx] = ind2sub( size(imLocalMax), seedInd );                
                theta = 0:0.1:(2*pi+0.1);
                
                cx = cx(:,ones(size(theta)));
                cy = cy(:,ones(size(theta)));
                rad = (sigmaValues(pixelScaleMap(seedInd)))';
                rad = rad(:,ones(size(theta)));
                
                theta = theta(ones(size(cx,1),1),:);                
                X = cx + cos(theta).* rad;
                Y = cy + sin(theta).* rad;

                imseriesmaskshow( imMultiscaleLoBGResponse, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Multiscale LoBG Filter' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
                
                figure, histogram( pixelScaleMap(seedInd) );
                title( 'Histogram of the cell scales (diameter) found in the image' );
                
            end
            
        case 3
            
            imLocalMax = locmax3d(imMultiscaleLoBGResponse, MaximaSuppressionSize);           
            imLocalMax = double(imLocalMax > 0);
            
            if flagDebugMode
                imseriesmaskshow( imMultiscaleLoBGResponse, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of Multiscale LoBG Filter' );
            end
            
    end

    % detect local intensity maxima as cell seed points
    imBlobLocations = imLocalMax;   
    if nargout > 1
        varargout{1} = imMultiscaleLoBGResponse;
    end     
    
end