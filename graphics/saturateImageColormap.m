function saturateImageColormap(handle,satPct)
%SATURATEIMAGECOLORMAP adjusts the colormap for the specified image so that it is saturated by the specified amount 
%
% success = saturateImageColormap(handle,satFrac)
%
%
% Hunter Elliott
% 8/24/2012
%

if nargin < 1 || isempty(handle)
    handle = gca;
end

if nargin < 2 || isempty(satPct)
    satPct = 1;
elseif satPct < 1
    %If they entered as a fraction.
    satPct = satPct * 100;
end

imDat = double(getimage(handle));
caxis(prctile(imDat(:),[satPct/2 (100-satPct/2)]));




