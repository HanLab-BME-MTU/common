function [maxXY,maxZY,maxZX,three]=computeMIPs(vol,ZXRatio,minInt,maxInt,varargin)
ip = inputParser;
ip.addParameter('stripeSize',8,@isnumeric);
ip.addParameter('tracks',[],@(t) isa(t,'Tracks'));
ip.addParameter('frameIdx',[],@isnumeric);
ip.parse(varargin{:});
p=ip.Results;

% set other parameters
stripeSize = p.stripeSize; % the width of the stripes in the image that combines all three maximum intensity projections
stripeColor = 0; %the stripe color, a number between 0 (black) and 1 (white).  (If you're not using all of the bit depth then white might be much lower than 1, and black might be higher than 0.)

ScaledZ=ceil(size(vol,3)*ZXRatio);
% find the maximum intensity projections
maxXY = (max(vol, [], 3));
maxZY = imresize((squeeze(max(vol, [], 2))),[size(vol,1) ScaledZ]);
maxZX = imresize((squeeze(max(vol, [], 1))),[size(vol,2) ScaledZ]);

% generate a single image with all three projections
threeTop = [maxXY, stripeColor*ones(size(vol,1), stripeSize), maxZY];
threeBottom = [maxZX', stripeColor*ones(ScaledZ, ScaledZ+stripeSize)];
three = [threeTop; stripeColor*ones(stripeSize, size(vol,2)+ScaledZ+stripeSize); threeBottom];

maxXY = uint8((2^8-1)*mat2gray(maxXY,double([minInt,maxInt])));
maxZY = uint8((2^8-1)*mat2gray(maxZY,double([minInt,maxInt])));
maxZX = uint8((2^8-1)*mat2gray(maxZX,double([minInt,maxInt])));
three = uint8((2^8-1)*mat2gray(three,double([minInt,maxInt])));
