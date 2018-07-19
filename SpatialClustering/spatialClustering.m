function [kr,lr,pcr]= spatialClustering(mpm,imsize,dist,plotMask,closureRadius,dilationRadius,doFill,sparse);
% spatialClustering puts together all functions required to calculate
% Ripley's functions with border correction
%           
%   INPUTS:
%           mpm:       mpm file containing (x,y) coordinates of points in 
%                      the image in succesive columns for different time 
%                      points
%           imsiz:     x,y-size of the image
%           dist =      bins for which correlation is calculated (e.g. 1:20
%           calculates the correlation at 1 pixel, 2pixels,...20 pixels)
%           plotMask =  1 to plot cell area mask and detections 0 to not
%           closureRadius =   radius for disk used in closure
%           dilationRadius = radius for dialation
%           doFill = 1 to fill mask, 0 to not fill
%           sparse: true to measure distances using graph based algorithm
%                   that creates sparse matrix (might be faster for large 
%                   data sets) 
%   OUTPUTS                    
%           kr:     Ripley's K-function 
%           lr:     Besag's L-function = sqrt(K(r))-r
%           pcr:    pair-correlation function
%           for all, every column contains the kr/lr function for one
%           frame of the mpm-file, the row values correspond to the
%           specified distances
%                                                
%

%
% D Nunez March 2, 2011


%interpret inputs
if nargin < 3 || isempty(dist)
    dist = 1:20;
end
if nargin < 4 || isempty(plotMask)
    plotMask = 0;
end
if nargin < 5 || isempty(closureRadius)
    closureRadius = 20;
end
if nargin < 6 || isempty(dilationRadius)
    dilationRadius = 5;
end
if nargin < 7 || isempty(doFill)
    doFill = 0;
end
if nargin < 8 || isempty(doFill)
    sparse = 0;
end

imsizS = [imsize(2) imsize(1)];

[areamask] = makeCellMaskDetections([mpm(:,1),mpm(:,2)],...
    closureRadius,dilationRadius,doFill,imsize,plotMask,[]);
normArea = bwarea(areamask);

corrFacMat = makeCorrFactorMatrix(imsizS, dist, 10, areamask');


[kr,lr,pcr]=RipleysKfunction(mpm,mpm,imsize,dist,corrFacMat,normArea);