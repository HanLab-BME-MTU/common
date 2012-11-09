function [imCellSeedPoints] = detect_cell_seeds_LoG( im, meanCellDiameter, varargin )
% Detects cell seed points as local maxima of the LoG filtered image
%
% The input intensity image is first filtered with a Laplacian of Gaussian (LoG)
% Filter and then the seed points are detected as local maxima in this
% filtered image.
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
    p.addRequired( 'meanCellDiameter', @(x) isscalar(x) );
    p.parse( im, meanCellDiameter );
    p.addParamValue( 'spacing', ones( 1, ndims(im) ), @(x) (isnumeric(x) && numel(x) == ndims(im)) );
    p.addParamValue( 'debugMode', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( im, meanCellDiameter, varargin{:} ); 
    
    spacing = p.Results.spacing;
    flagDebugMode = p.Results.debugMode;
    meanCellDiameterImsp = meanCellDiameter ./ spacing;

    % Apply Laplacian of Gaussian (LoG) filter
    sigma = 0.5 * meanCellDiameter / sqrt( ndims(im) );
    imLoG = -1 * filterLoGND( im, sigma, 'spacing', spacing, 'UseNormalizedDerivatives', true );

    % locate local intensity maxima in gaussian blurred image
    MaximaSuppressionSize = round( meanCellDiameterImsp );
    evenind = (mod( MaximaSuppressionSize, 2 ) == 0);
    MaximaSuppressionSize( evenind ) = MaximaSuppressionSize( evenind ) + 1;    
    
    switch ndims( im ) 
        
        case 2 
            
            imLocalMax = locmax2d(imLoG, MaximaSuppressionSize, 1);            
    
            if flagDebugMode

                [cy, cx] = find( imLocalMax > 0 );                
                theta = 0:0.1:(2*pi+0.1);
                cx = cx(:,ones(size(theta)));
                cy = cy(:,ones(size(theta)));
                cellRadius = (0.5 * meanCellDiameter);
                rad = cellRadius * ones( size(cx) );
                theta = theta(ones(size(cx,1),1),:);                
                X = cx + cos(theta).* rad;
                Y = cy + sin(theta).* rad;
                
                imseriesmaskshow( imLoG, imdilate( imLocalMax, strel('disk',3) ) ); 
                set( gcf, 'Name', 'Local Maxima Overlayed on the response of LoG Filter' );
                hold on;               
                    line(X', Y', 'Color', 'g', 'LineWidth', 1.5);                
                hold off;
                
            end
            
        case 3
            
            imLocalMax = locmax3d(imLoG, MaximaSuppressionSize);           

            if flagDebugMode
                imseriesmaskshow( imLoG, imdilate(imLocalMax, ones(3,3,3)) ); 
                set( gcf, 'Name', 'Local Maxima overlayed on the response of LoG Filter' );
            end
            
    end

    % detect local intensity maxima as cell seed points
    imCellSeedPoints = imLocalMax;

end