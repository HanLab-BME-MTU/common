function [ h ] = imshowpair_scrollfade( varargin )
%IMSHOWPAIR_SCROLLFADE Use mouse scroll wheel to fade out a channel
%
% INPUT
% imshowpair_scrollfade( imshowpair_parameters );
% imshowpair_scrollfade( {imshowpair_parameters}, scroll_parameters);
% imshowpair_scrollfade( rgb, scroll_parameters);
% 
% rgb is an image with size(rgb,3) == 3 such as created by imfuse
%
% PARAMETERS (for scrolling)
% ScrollSpeed - Amount that each mouse wheel click increases or decreases
%     the intensity. Use negative values to invert the direction of wheel.
%     default: 20
% ScrollMax - How many times the original maximum to enhance the channel
%     default: 1
% Channels - Which channels to fade in and out
%     default: 2
%
% OUTPUT
% h - handle to matlab.graphics.primitive.Image by imshow/imshowpair
%
% Mark Kittisopikul, February 2016
% Jaqaman Lab
% UT Southwestern

params = {};

rgb = [];

if(iscell(varargin{1}))
    passthrough = varargin{1};
    if(nargin > 1)
        params = varargin(2:end);
    end
elseif(size(varargin{1},3) == 3)
    rgb = varargin{1};
    params = varargin(2:end);
else
    passthrough = varargin;
end

ip = inputParser;
ip.addParameter('ScrollSpeed',20,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('ScrollMax',1,@(x) isscalar(x) && isnumeric(x));
ip.addParameter('Channels',2,@isnumeric);
ip.parse(params{:});

fh = figure;
if(~isempty(rgb))
    h = imshow(rgb,[]);
else
    h = imshowpair(passthrough{:});
    rgb = double(h.CData);
end

fh.WindowScrollWheelFcn = @scrollFade;

orgMax = max(joinColumns(rgb(:,:,ip.Results.Channels)));
intensity_max = orgMax;


    function scrollFade(~,cb)
        intensity_max = intensity_max + cb.VerticalScrollCount*ip.Results.ScrollSpeed;
        intensity_max = max(intensity_max,0);
        intensity_max = min(intensity_max,ip.Results.ScrollMax*orgMax);
        scaling = ones(1,3);
        scaling(ip.Results.Channels) = intensity_max/orgMax;
        I = uint8(bsxfun(@times,rgb,shiftdim(scaling,-1)));
        h.CData = I;
    end

end

