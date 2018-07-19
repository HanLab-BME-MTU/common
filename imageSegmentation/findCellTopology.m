function locmax = findCellTopology(image,scale,varargin)
% uses gaussian filter and finds local maxima in order to localize cells in
% an image
%
%   Inputs: image - input image (required)
%           scale - std of gaussian as well as neighborhood for local
%                   maxima in pixels (default: 50)
%           plot - true to plot image and result side by side
%                   (default:false)
%           padFactor - image is padded by this factor multiplied by the 
%                   scalenumber of pixels (default: 1)
%           reSizeFactor - image is downsampled to a square of side length
%                   equal to this factor
%
%   Outputs: locmax - mask where center of cells are marked by a 1
%
%
% dnunez 2012 October 26

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('image', @isnumeric);
ip.addParamValue('scale', 50 ,@isscalar);
ip.addParamValue('plot', false, @islogical);
ip.addParamValue('padFactor', 1, @isinteger);
ip.addParamValue('reSizeFactor', 100, @isinteger);
ip.parse(image,scale,varargin{:});
scale = ip.Results.scale;
padFactor = ip.Results.padFactor;
reSizeFactor = ip.Results.reSizeFactor;

%low pass filter image
filteredImage = filterGauss2D(image,scale);

%pad symmetrically by padFactor*scale pixels on all sides 
filteredImage = padarray(filteredImage,[round(padFactor*scale) round(padFactor*scale)],'symmetric','both');

%downasmple filtered image to reSizeFactorXreSizeFactor
filteredImageResized = resizem(filteredImage,[reSizeFactor reSizeFactor]);
%find local maxima
locmax = locmax2d(filteredImageResized,round(scale*reSizeFactor/size(image,1)),1);
%upsample mask with local maima
locmax = resizem(locmax,size(filteredImage));

%remove padding
filteredImage = filteredImage(padFactor*scale+1:end-padFactor*scale,padFactor*scale+1:end-padFactor*scale);
locmax = locmax(padFactor*scale+1:end-padFactor*scale,padFactor*scale+1:end-padFactor*scale);

%plot results
if ip.Results.plot
    
    figure
    subplot(1,2,1), imagesc(image)
    subplot(1,2,2), imagesc(image), hold on
    se = strel('rectangle',[round(scale/2) round(scale/2)]);
    locmaxd = imdilate(locmax,se);
    boundaries = bwboundaries(locmaxd);
    for i = 1:length(boundaries)
        plot(boundaries{i}(:,2),boundaries{i}(:,1),'r')
    end
    
end
