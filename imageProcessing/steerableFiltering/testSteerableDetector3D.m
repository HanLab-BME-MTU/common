% 
% The following script tests the steerableDetector3D on a synthetic 3D volume 
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

    numLines = 100;
    
    meanLineWidth = 2;
    stdLineWidth = 0;

    sigmaSteerableDetector = meanLineWidth;
    
%*************************************************************************

im = zeros( 200, 200, 30 );

% normal distributons for foreground and background pixels
fgMeanVar = [ 200, 20 ];
bgMeanVar = [ 180, 20 ];
fgGmObj = gmdistribution( fgMeanVar(:,1), reshape(fgMeanVar(:,2), [1,1,size(fgMeanVar,1)]) );
bgGmObj = gmdistribution( bgMeanVar(:,1), reshape(bgMeanVar(:,2), [1,1,size(bgMeanVar,1)]) );
im(:) = random(bgGmObj,numel(im));

% generate lines randomly
fprintf( '\nGenerating %d Lines in Random Orientations: \n', numLines );

imsize = size( im );
[X,Y,Z] = meshgrid(1:imsize(2), 1:imsize(1), 1:imsize(3));

try
    ptIm = gpuArray( [ X(:), Y(:), Z(:) ] );
catch err
    ptIm = [ X(:), Y(:), Z(:) ]; % you probably dont have a gpu, the code will be slow
end

lineMask = false( imsize );

for i = 1:numLines   
    
    fprintf( ' %.3d ', i );
    
    if mod( i, 15 ) == 0
        fprintf( '\n' );
    end
    
    % generate a random point in volume
    ptRandInd = floor( rand * numel(im) );
    ptRef = ptIm( ptRandInd, : );
    
    % generate random orientation vector
    thetax = rand * pi;
    thetaxy = (0.05 + 0.05 * randn ) * pi;
    v = [ cos(thetaxy) * cos(thetax), cos(thetaxy) * sin(thetax), sin(thetaxy) ];
    
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

% display
imseriesshow( im );
set( gcf, 'Name', 'Volume with Randomly Generated Lines' );

[res, theta, nms] = steerableDetector3D(im, 1, sigmaSteerableDetector);

imseriesshow( res );
set( gcf, 'Name', 'Response of Steerable Detector' );

imseriesshow( nms );
set( gcf, 'Name', 'Result of Non-maximal suppression' );
