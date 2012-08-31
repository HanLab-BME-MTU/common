% 
% The following script tests the steerableDetector on a synthetic 2D image
% with lines drawn in random orientations
% 
% Author: Deepak Roy Chittajallu (Created on Aug 31, 2012)
%

clc
clear
close all

%*************************************************************************
%                               PARAMETERS
%*************************************************************************

    numLines = 20;
    
    meanLineWidth = 2;
    stdLineWidth = 0.0;
    
    sigmaSteerableDetector = meanLineWidth;
    orderSteerableDetector = 4;

%*************************************************************************

im = zeros( 512, 512 );

% normal distributons for foreground and background pixels
fgMeanVar = [ 200, 20 ];
bgMeanVar = [ 180, 20 ];
fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );
im(:) = random(bgGmObj,numel(im));

% generate lines randomly
fprintf( '\nGenerating %d Lines in Random Orientations: \n', numLines );

imsize = size( im );
[X,Y] = meshgrid(1:imsize(2), 1:imsize(1));

try
    ptIm = gpuArray( [ X(:), Y(:) ] );
catch err
    ptIm = [ X(:), Y(:) ]; % you probably dont have a gpu, the code will be slow
end

lineMask = false( imsize );

for i = 1:numLines   
    
    fprintf( ' %.3d ', i );
    
    if mod( i, 15 ) == 0
        fprintf( '\n' );
    end
    
    % generate a random point in the image
    ptRandInd = floor( rand * numel(im) );
    ptRef = ptIm( ptRandInd, : );
    
    % generate random orientation vector
    thetax = rand * pi;
    v = [ cos(thetax), sin(thetax) ];
    
    % generate random line width
    randLineWidth = meanLineWidth + stdLineWidth * rand;
    
    if randLineWidth < 0 
        randLineWidth = meanLineWidth;
    end
    
    % generate line
    refvec = ptIm - repmat( ptRef, numel(im), 1 );
    perp = refvec - (refvec * v') * v;
    sqDist = sum( perp .* perp, 2 );
    lineMask( sqDist <= randLineWidth^2 ) = true;
    
end

fprintf( '\n' );

im( lineMask ) = random(fgGmObj, numel( find( lineMask ) ));

% Run steerable detector to enhance the lines
fprintf( '\nRunning steerable detector to enhance curves on %d x %d sized image ...\n', imsize(2), imsize(1) );
tic
[res, theta, nms, rotations] = steerableDetector(im, orderSteerableDetector, sigmaSteerableDetector);
timeElapsed = toc;
fprintf( '\n\tIt took %.2f seconds\n', timeElapsed );

% display
imseriesshow( im );
set( gcf, 'Name', 'Image with Randomly Generated Lines' );

imseriesshow( res );
set( gcf, 'Name', 'Response of Steerable Detector' );

imseriesshow( nms );
set( gcf, 'Name', 'Result of Non-maximal suppression' );
