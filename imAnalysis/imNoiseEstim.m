function [nse,thresh,grad]=imNoiseEstim(img,c)
%IMNOISEESTIM returns an estimate of the image noise
% The algorithms is described in Voorhees and Poggio, 
%                                Darpa Image Understanding Workshop 892-899, 1987.
%
% SYNOPSIS [nse,thresh,grad]=imNoiseEstim(img,c)
%
% INPUT img : greyvalue image (uint8 or double accepted)
%       c   : confidence probability on which the threshold is set (optional)
%             default 0.99
%
% OUTPUT nse    : noise estimate (in grayvalues)
%        thresh : threshold for any type of low evidence suppression
%        grad   : buffer with the gradients (allows one to display the 
%                 approximately Rayleigh distributed gradient field)
%
% Matlab-function by KQ 2003 (replaces C-version by GD)


if(nargin == 1)
   c = 0.01;
end;

if(~(isa(img,'double') | isa(img,'uint8')))
   error('Invalid data type entered');
end;

if(isa(img,'uint8'))
   auxImg = double(img)/255;
   [nse,thresh,grad]=NoiseEstim(auxImg,c);
   nse = 255*nse;
   thresh = 255*nse;
   grad = uint8(round(grad*255));
else
   [nse,thresh,grad]=NoiseEstim(img,1-c);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTION NoiseEstim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nse,thresh,grad]=NoiseEstim(auxImg,c)

% Convolve the image with the central difference kernel to 
% obtain the image gradient; calculate gradient magnitude    
diff=[-0.5;0;0.5];               
ax = conv2(auxImg, diff', 'same');  
ay = conv2(auxImg, diff, 'same'); 
mag = sqrt((ax.*ax) + (ay.*ay)); 

% sort gradient magnitudes in ascending order, 
% reduce to non-zero elements
mag = sort(mag(:));
magr = mag(min(find(mag)):end);

% Gradient magnitude is (approximately) Raleigh-distributed. Look 
% for first mode of the distribution in cumulative histogram
winSize = floor([.1:.05:.5]*length(magr));
diffMat = [1:floor(length(magr)/2)]';
for i = 1:length(winSize)
    diffMat = [diffMat, diffMat(:,1)+winSize(i)]; % 10 cols
end
data = magr(diffMat);% 10 cols
diffData = []; 
for i = 2:length(winSize)+1
    diffData = [diffData, data(:,i)-data(:,1)]; % 9 cols
end
idx = min(diffData); %???
maxlist = [];
for i = 1:length(winSize)
    row = find(diffData(:,i)==idx(i));
    max = (data(row,i+1)+data(row,1))/2;
    maxlist = [maxlist;max];
    % middle of the window position where difference between both borders is minimal
end
mode = mean(maxlist);

% according to the paper: noise = mode / sqrt(2);
% BUT: since the gradient has been computed over the distance 2 instead of 1
% (lowpass filtering) which is gives a "noise gradient" of a factor 2 too 
% small (the probability to have intensity variation between two locations in an 
% unstructured, noisy image is independent of the distance between the 2 samples
% used to compute the gradient). thus, instead of division by sqrt(2) one has 
% to multiply the mode with 2/sqrt(2)=sqrt(2).
	
nse = sqrt(2)*mode;
thresh = sqrt(-2.0*log(c))*mode;
grad = mag;



