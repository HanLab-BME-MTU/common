function saturateImageColormap(handle,satPct)
%SATURATEIMAGECOLORMAP adjusts the colormap for the specified image so that it is saturated by the specified amount 
%
% success = saturateImageColormap(handle,satPct)
%  
%   handle - axes handle with image to saturate
%   satPct - the percentage of pixels to saturate. (Half this amount will
%            be saturated on the low end, half on the high end)
%
% Hunter Elliott
% 8/24/2012
%

if nargin < 1 || isempty(handle)
    handle = gca;
end

if nargin < 2 || isempty(satPct)
    satPct = 1;
end

imDat = double(getimage(handle));
caxis(prctile(imDat(:),[satPct/2 (100-satPct/2)]));




