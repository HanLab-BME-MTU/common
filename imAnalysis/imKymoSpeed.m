function v = imKymoSpeed(kym,bw)
%imKymoSpeed: Automatically detect in a kymograph the average line slope that 
%             gives the speed of a flow.
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
      'the band width, ''bw''. See help imKymoSpeed.'])
end

if numFrames < 2
   error('There have to be at least two frames in the kymograph.');
end

%A threathold is set to determine if the kymograph is distinguishable to get
% the speed. Basically, it theck if the difference in the scores for all
% speeds is big enough so that the calculated speed is acceptable. If the
% difference is too small, 'NaN' is returned. This value is set empirically.
acceptThreshold = 0.02;

%The length of the smallest line segment the program attempts to detect.
% The length the line is in terms of the number of frames it spans.
lineSegLen = 3;

%The length of the line (or curve) in the movie image. It is also the
% horizontal length of the kymograph.
kymHLen = size(kym,2);

%We choose the maximum speed to be 'kymHLen-1'. At this speed, it means 
% that the movie image will flow the whole length of the selected curve for
% kymograph analysis.
%speedLimit  = kymHLen-1;
speedLimit  = ceil(kymHLen/2);

%Reshape 'kym' so that we can correlate between frames.
bandKym = squeeze(sum(reshape(kym.',kymHLen,bw,numFrames),2)).';

%Scale the intensity of the kymograph to be from 0.1 to 1 to increase contrast
% and to increase the difference between scores.
minKymI = min(bandKym(:));
maxKymI = max(bandKym(:));

if maxKymI-minKymI ~= 0
   bandKym = (bandKym-minKymI)/(maxKymI-minKymI)*1 + 0.1;
end

%The square of the intensity matrix for the use of normalization later.
bandKymP2 = bandKym.*bandKym;

%Calculate a score for each discretized speed or line slope in the range
% (-speedLimit,speedLimit). The step size if one pixel.
% v     : flow speed or the slope of lines in the kymograph.
% score : Correlation score for each line slope.
score  = zeros(2*speedLimit+1,1);
%Negative speed.
for v = -speedLimit:-1
   k = v+speedLimit+1;

   %We calculate scores for all line segments of length 'lineSegLen' and in
   % the direction of slope, v.
   %First, check all the lines that start from each pixel on the top edge of 
   % the banded kymograph, 'bandKym'.
   numLineSegs = 0;
   for j = 2:kymHLen
      %The total length of the line in terms of frames.
      totalLen = min(numFrames,floor((j-1)/-v));
      for f2 = lineSegLen:totalLen
         numLineSegs = numLineSegs + 1;
         f1 = f2-lineSegLen+1;
         %The frames over which the line segment spans.
         f  = f1:f2;

         %The column indices of the line segment between frame 'f1' and 'f2'.
         col = j+(f-1)*v;

         %The indices of the line segment in the banded kymograph matrix,
         %'bandKym' when it is viewed as a 1D array.
         ind = (col-1)*numFrames+f;

         %The product of the elements on the line segment adds to the total
         % score for the current slope, v.
         score(k) = score(k) + prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))));
      end
   end

   %Then, check all the lines that start from each pixel on the right edge of 
   % the banded kymograph, 'bandKym'.
   for j = 2:numFrames
      totalLen = min(numFrames-j+1,floor(kymHLen/-v));
      for f2 = lineSegLen:totalLen
         numLineSegs = numLineSegs+1;
         f1 = f2-lineSegLen+1;
         f  = f1:f2;

         col = kymHLen+(f-1)*v;
         ind = (col-1)*numFrames+f+j-1;

         score(k) = score(k) + prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))));
      end
   end

   if numLineSegs ~= 0
      score(k) = score(k)/numLineSegs;
   end

   %Each element of the integrand.
   %corrM = bandKym(1-v:end,:,1:numFrames-1).*bandKym(1:end+v,:,2:numFrames);

   %Calculate the norm of the intensity of the overlapped band after shift.
   % The shift is '-v'.
   %normIL(:,1-v) = sqrt(squeeze(sum(sum(bandKymP2(1-v:end,:,:),2),1)));
   %normIR(:,1-v) = sqrt(squeeze(sum(sum(bandKymP2(1:end+v,:,:),2),1)));

   %Normalized correlation.
   %corrI    = squeeze(sum(sum(corrM,2),1))./normIL(1:numFrames-1,1-v)./ ...
   %           normIR(2:numFrames,1-v);
   %corrI    = squeeze(sum(sum(corrM,2),1))/(kymHLen+v);
   %score(k) = sum(corrI);
end

%Zero Speed: we only have lines that start from pixels on the top edge of
% 'bandKym' and the total length is always 'numFrames'.
k = speedLimit+1;
numLineSegs = 0;
for j = 1:kymHLen
   for f2 = lineSegLen:numFrames
      numLineSegs = numLineSegs+1;
      f1 = f2-lineSegLen+1;
      f  = f1:f2;

      score(k) = score(k) + prod(bandKym(f,j)/sqrt(sum(bandKymP2(f,j))));
   end
end

if numLineSegs ~= 0
   score(k) = score(k)/numLineSegs;
end

%Positive speed.
for v = 1:speedLimit
   k = v+speedLimit+1;

   numLineSegs = 0;
   %First, check all the lines that start from each pixel on the top edge of 
   % the banded kymograph, 'bandKym'.
   for j = 2:kymHLen
      %The total length of the line in terms of frames.
      totalLen = min(numFrames,floor((j-1)/v));
      for f2 = lineSegLen:totalLen
         numLineSegs = numLineSegs + 1;
         f1 = f2-lineSegLen+1;
         %The frames over which the line segment spans.
         f  = f1:f2;

         %The column indices of the line segment between frame 'f1' and 'f2'.
         col = kymHLen-j+1+(f-1)*v;

         %The indices of the line segment in the banded kymograph matrix,
         %'bandKym' when it is viewed as a 1D array.
         ind = (col-1)*numFrames+f;

         %The product of the elements on the line segment adds to the total
         % score for the current slope, v.
         score(k) = score(k) + prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))));
      end
   end

   %Then, check all the lines that start from each pixel on the left edge of 
   % the banded kymograph, 'bandKym'.
   for j = 2:numFrames
      totalLen = min(numFrames-j+1,floor(kymHLen/v));
      for f2 = lineSegLen:totalLen
         numLineSegs = numLineSegs+1;
         f1 = f2-lineSegLen+1;
         f  = f1:f2;

         col = 1+(f-1)*v;
         ind = (col-1)*numFrames+f+j-1;

         score(k) = score(k) + prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))));
      end
   end

   if numLineSegs ~= 0
      score(k) = score(k)/numLineSegs;
   end

   %corrM = bandKym(1:end-v,:,1:numFrames-1).*bandKym(v+1:end,:,2:numFrames);

   %Normalized correlation.
   %corrI    = squeeze(sum(sum(corrM,2),1))./normIR(1:numFrames-1,1+v)./ ...
   %           normIL(2:numFrames,1+v);
   %corrI    = squeeze(sum(sum(corrM,2),1))/(kymHLen-v);
   %score(k) = sum(corrI);
end

[maxScore,k] = max(score);
v = k-speedLimit-1;

%Check if the difference in the score is significant enough to accept the
% identified speed.
if maxScore-sum(score)/length(score) < maxScore*acceptThreshold
   v = NaN;
   return;
end

%Fine tune: subpixel velocity tracking. We first interpolate the speed at
% sampling pixels with cubic spline and then find the maximum.
vv = [1:length(score)]-speedLimit-1;
pp = csape(vv,score);
v  = fminbnd(@scoreFun,max(vv(1),v-1),min(vv(end),v+1),[],pp);

function f = scoreFun(v,pp)

f = -fnval(pp,v);
