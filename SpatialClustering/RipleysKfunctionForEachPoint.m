function [kr,lr,pcr,glr]=RipleysKfunctionForEachPoint(mpm,imsize,varargin)
% RipleysKfunction calculates Ripley's K-function for a given MPM,
% allowing cross-corrlation between two MPMs
% SYNOPSIS  [kr,lr,pcr]=RipleysKfunction(mpm,imsiz,dist,corrFacMat, normArea);
%
%INPUT  mpm                 : matrix of point coordinates [x1 y1; x2 y2;....xn yn]
%       imsize              : image size of the form [numVerticalPixels numHorizontalPixels]
%       dist(optional)      : vector of distances at which Ripley's K function
%                             (i.e. number of neighbors )is sampled (1:20 is default)
%       mask(optional)      : user provided binary mask used to calculate area
%                             and to correct point density
%       calculateMask(opt)  : true to calculate mask from mpm positions, false
%                             otherwise
%       closureRadius(opt)  : radius for disk used in closure for calcuating
%                             mask
%       dilationRadius(opt) : radius for dialation in calculating mask
%       doFill(optional)    : true to fill calculates mask, false to not fill
%       plotMask(optional)  : true to display calculated mask, false to not
%
%OUTPUT
%           kr:     Ripley's K-function
%           lr:     Besag's L-function = sqrt(K(r))-r
%           pcr:    pair-correlation function
%           glr:     Getis' L-function = sqrt(K(r))
%           for all, every column contains the kr/lr function for one
%           point of the mpm-file, the row values correspond to the
%           specified distances
%
%
%Daniel Nunez, April 2012


ip = inputParser;
ip.CaseSensitive = false;
%simulation params
ip.addRequired('mpm', @isnumeric);
ip.addRequired('imsize', @isnumeric);
ip.addParamValue('dist', 1:20, @isnumeric);
ip.addParamValue('mask', [], @islogical);
ip.addParamValue('calculateMask', false, @islogical);
ip.addParamValue('closureRadius', 20, @islogical);
ip.addParamValue('dilationRadius', 5, @islogical);
ip.addParamValue('doFill', false, @islogical);
ip.addParamValue('plotMask', false, @islogical);
ip.parse(mpm, imsize, varargin{:});
dist = ip.Results.dist;
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
lr = nan(length(dist),size(mpm,1));
kr = lr;
pcr = lr;
glr = lr;
for i = 1:length(mpm)
    [kr(:,i),lr(:,i),pcr(:,i),glr(:,i)]=RipleysKfunction(mpm([1:i-1 i+1:length(mpm)],:),mpm(i,:),imsizS,dist,corrFacMat,normArea);
end