function outputImage = multicolorImage(images,cmap,norm)
%MULTICOLORIMAGE creates an overlay color image out of multiple 2D images
%
% SYNOPSIS: outputImage = multicolorImage(images,colormap)
%
% INPUT images: xPix-by-yPix-by-nImages stack of images
%		cmap: (opt) nImages-by-3 colormap (see help colormap), or string of
%		      color codes ('rgb' resolves to red, green, blue)
%             Default: Jet
%       norm: (opt) If 1, individual images will be normed to 0...1.
%             Default: 1
%
% OUTPUT outputImage: multicolor overlay image
%
% REMARKS:
%
% created with MATLAB ver.: 7.5.0.342 (R2007b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 10-Mar-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=========================
%% TEST INPUT
%=========================

% check images
if nargin == 0 || isempty(images)
    error('multicolorImage needs at least one nonempty input argument')
end
if ndims(images) > 3
    error('muticolorImage can only handle stacks of 2d images')
end
imSize = size(images);
if length(imSize) < 3
    nImages = 1;
else
    nImages = imSize(3);
end
imSize = imSize(1:2);

% convert to double
images = double(images);


% check for colormap
if nargin < 2 || isempty(cmap)
    % there are better maps than jet, but leave it at this for the moment
    cmap = colormap(jet(nImages));
end
% resolve rgb, rgbw etc
if ischar(cmap)
    if length(cmap) == nImages
        tmp = zeros(nImages,3);
        for i=1:nImages
            tmp(i,:) = colorCode2rgb(cmap(i));
        end
        cmap = tmp;
    else
        error('colormap needs as many colors as there are images')
    end
end
if ~all(size(cmap) == [nImages,3])
    error('colormap must be nImages-by-3 array')
end

if nargin < 3 || isempty(norm)
    norm = true;
end

%========================


%========================
%% CREATE IMAGE
%========================

% create output image
outputImage = zeros([imSize,3]);

% add individual images to end result. A lot of overlay may lead to
% "saturation". I'm not correcting for that for now.
for i=1:nImages
    if norm
        tmp = norm01(images(:,:,i));
    else
        tmp = images(:,:,i);
    end
    outputImage = outputImage + ...
        cat(3,tmp*cmap(i,1),...
        tmp*cmap(i,2),...
        tmp*cmap(i,3));
end

if nargout == 0
    figure,imshow(outputImage);
end