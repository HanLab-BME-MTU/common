function [C, H] = plotLocalNeighborhoodDensityMap(mpm,imsize,varargin)

%makes 'Getis Map' of point processes to identify clusters
%SYNOPSIS [C, H] = plotLocalNeighborhoodDensityMap(mpm,imsize,varargin))
%
%INPUT  mpm             : matrix of point coordinates [x1 y1; x2 y2;....xn yn]
%       imsize          : image size of the form [numVerticalPixels numHorizontalPixels]
%       dist(optional)  : vector of distances at which Ripley's K function
%                         (i.e. number of neighbors )is sampled (1:20 is default)
%       scale(optional) : intensity of clustering to map
%       lr(optional)    : can specify Besag's L-function to use as
%                         lr(length(dist),numOfPoints), use
%                         RipleysKfunction.m
%       mask(optional)      : user provided binary mask used to calculate area
%                             and to correct point density
%       calculateMask(opt)  : true to calculate mask from mpm positions, false
%                             otherwise
%       closureRadius(opt)  : radius for disk used in closure for calcuating
%                             mask
%       dilationRadius(opt) : radius for dialation in calculating mask
%       doFill(optional)    : true to fill calculates mask, false to not fill
%       plotMask(optional)  : true to display calculated mask, false to not
%       function(optional)  : specify which function to use kr, lr, gcr, or
%                             glr (see RipleysKfunction.m)
%       sparse          : true to measure distances using graph based algorithm
%                         that creates sparse matrix (might be faster for large 
%                         data sets) 
%
%OUTPUT C               : contour matrix C as described in CONTOURC
%       H               : handle H to a contourgroup object
%
%Daniel Nunez, April 2012
%
%Note: This algorithm essentially assigns the value of lr evaluated at a
%specified distance (the scale) for each possition specified in the mpm.
%For sparce data, additional random points are created for which the lr is
%calculated. These additional points do not factor into the lr measured
%for other points, and therefore serve as a way to sample the neighbor
%density localy. This was addapted from Getis and Franklin (1987).


ip = inputParser;
ip.CaseSensitive = false;
%simulation params
ip.addRequired('mpm', @isnumeric);
ip.addRequired('imsize', @isnumeric);
ip.addParamValue('dist', 1:20, @isnumeric);
ip.addParamValue('scale', 1, @isscalar);
ip.addParamValue('lr', [], @isnumeric);
ip.addParamValue('mask', [], @islogical);
ip.addParamValue('calculateMask', false, @islogical);
ip.addParamValue('closureRadius', 20, @islogical);
ip.addParamValue('dilationRadius', 5, @islogical);
ip.addParamValue('doFill', false, @islogical);
ip.addParamValue('plotMask', false, @islogical);
ip.addParamValue('function', 'lr', @(var)any(strcmp(var,{'kr' 'lr' 'gcr' 'glr'})));
ip.addParamValue('sparse', false, @islogical);
ip.parse(mpm, imsize, varargin{:});
dist = ip.Results.dist;
lr = ip.Results.lr;
areamask  = ip.Results.mask;
imsizS = [imsize(2) imsize(1)];


%MAKE MASK
if ip.Results.calculateMask
    [areamask] = makeCellMaskDetections(mpm,ip.Results.closureRadius,...
        ip.Results.dilationRadius,ip.Results.doFill,imsize,ip.Results.plotMask,[]);
end

if ~isempty(areamask)
    %CALCULATE NORMALIZED AREA FROM MASK
    normArea = bwarea(areamask);
    % CREATE CORRECTION FACTOR MATRIX FOR THIS MOVIE using all objects
    corrFacMat  = makeCorrFactorMatrix(imsizS, dist, 10, areamask');
else
    corrFacMat = [];
    normArea = [];
end


%Calculate L-ripley
if isempty(lr)
    lr = nan(length(dist),size(mpm,1));
    for i = 1:length(mpm)
        [k,l,g,gl]=RipleysKfunction(mpm([1:i-1 i+1:length(mpm)],:),mpm(i,:),imsizS,dist,corrFacMat,normArea,ip.Results.sparse);
        
        if strcmp(ip.Results.function,'kr')
            lr(:,i) = k;
        elseif strcmp(ip.Results.function,'lr')
            lr(:,i) = l;
        elseif strcmp(ip.Results.function,'gcr')
            lr(:,i) = g;
        elseif strcmp(ip.Results.function,'glr')
            lr(:,i) = gl;
        end
        
    end
end

%make map out of L values for each point at a given scale (i.e., for each
%point assign it's L(d) for a given d)
M = zeros(imsize);
for i = 1:size(lr,2)
    M(max(1,round(mpm(i,2))),max(1,round(mpm(i,1)))) = lr(ip.Results.scale,i);
end


%add random points to make sure plot fills out nicely
percentPoints = 4;
mpmRandSample = repmat(imsize,percentPoints*imsize(1)*imsize(2),1).*...
    rand(percentPoints*imsize(1)*imsize(2),2);
lrRandSample = nan(length(dist),size(mpm,1));
for i = 1:length(mpmRandSample)
    [~,lrRandSample(:,i)]=RipleysKfunction(mpm,mpmRandSample(i,:),imsizS,dist,corrFacMat,normArea);
end
for i = 1:size(lrRandSample,2)
    M(max(1,round(mpmRandSample(i,2))),max(1,round(mpmRandSample(i,1)))) = lrRandSample(ip.Results.scale,i);
end


%plot contours
%V = linspace(min(M(:)),max(M(:)),10);
figure
[C H] = contourf(M);
clabel(C,H);
hold on
plot(mpm(:,1),mpm(:,2),'m.')