function v = imKymoSpeed(kym,bw)
%imKymoSpeed: Automatically detect in a kymograph the line whose slope gives 
%             the speed of a flow.
%
% SYNTAX: v = imKymoSpeed(kym,bw)
%    kym : The kymograph dissected and stacked from the movie of images. It is
%          a two dimentional image array whose horizontal size (in pixels) is
%          the length of the line selected in the movie image and whose
%          vertical size is 'bw*numFrames'.
%    bw  : The band width of the tube around a line in the movie images whose
%          kymograph is 'kym'.
%    v   : The speed: slope of the detected line in the kymograph.
%

%Check input.
numFrames = size(kym,1)/bw;
if ceil(numFrames) ~= numFrames
   error(['The size of the kymograph, ''kym'' does not match' ...
      'the band width, ''bw''. See help imKymoSpeed.']
end

if numFrames < 2
   error('There have to be at least two frames in the kymograph.');
end

%The length of the line (or curve) in the movie image.
kymHLen = size(kym,2);

%We choose the maximum speed to be 'kymHLen-1'. At this speed, it means 
% that the movie image will flow the whole length of the selected curve for
% kymograph analysis.
speedLimit  = kymHLen-1;

%Calculate the target function for each discretized speed in the range
% (-speedLimit,speedLimit). The step size if one pixel.
%Reshape 'kym' so that we can correlate between frames.
bandKym = reshape(Kym.',kymHLen,bw,numFrames);
%The total intensity of the band for each frame.
totalI = sum(sum(bandKym,2),1);
fCorr   = zeros(2*speedLimit-1,1);
for v = -speedLimit+1:0
   %Negative speed.
   k     = v+speedLimit;
   corrM = bandKym(k+1:end,:,1:numFrames-1).*bandKym(1:end-k,:,2:numFrames);

   %Normalized correlation.
   corrI    = sum(sum(corrM,2),1)./totalI(1:numFrames-1)./totalI(2:numFrames);
   fCorr(k) = sum(corrI);
end

for v = 1:speedLimit-1
   %Positive speed.
   k     = v+speedLimit;
   corrM = bandKym(1:end-v,:,1:numFrames-1).*bandKym(v+1:end,:,2:numFrames);

   %Normalized correlation.
   corrI    = sum(sum(corrM,2),1)./totalI(1:numFrames-1)./totalI(2:numFrames);
   fCorr(k) = sum(corrI);
end

[corrMax,k] = max(fCorr);
v = k-speedLimit;
