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

%Try 'imLineDetect'.
%opt.scales = 1:3;
%opt.linetype = 7;
%[binMask,resp,ori] = imLineExtract(double(kym),opt);

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
% the speed. Basically, it theck if the difference between the max score and
% the average is significantly bigger than the deviation of the scores for all
% speeds. For example, if acceptThreshold is 2, it means that the difference
% has to be 2 times bigger than the deviation. Otherwise, 0 is returned. 
% This value is set empirically.
acceptThreshold = 2;

%The length of the smallest line segment the program attempts to detect.
% The length the line is in terms of the number of frames it spans.
lineSegLen = 3;

lineThreshold = 0;

%The length of the line (or curve) in the movie image. It is also the
% horizontal length of the kymograph.
kymHLen = size(kym,2);

%Reshape 'kym' so that we can correlate between frames.
bandKym = squeeze(sum(reshape(kym.',kymHLen,bw,numFrames),2)).';

%Scale the intensity of the kymograph to be from 0.1 to 1 to increase contrast
% and to increase the difference between scores.
minKymI = min(bandKym(:));
maxKymI = max(bandKym(:));

if maxKymI-minKymI ~= 0
   bandKym = (bandKym-minKymI)/(maxKymI-minKymI)*1 + 0.1;
end

%Interpolate to half pixels so that we can sample half pixel velocity. After
% interpolation, we treat half pixels as integers. In the end, the slope has
% to be divided by 2 to get the real velocity.
ppKym = csape([1:kymHLen],bandKym);
bKymFine = fnval(ppKym,[1:0.5:kymHLen]);
bandKym = bKymFine;
kymHLen = size(bandKym,2);

%The square of the intensity matrix for the use of normalization later.
bandKymP2 = bandKym.*bandKym;

%We choose the maximum speed so that at this speed, a speckle will 
% flow 1/3 (if lineSeg = 3) the whole length of the selected curve for
% kymograph analysis. It is multiplied by 2 because we interpolate to half
% pixels.
%speedLimit  = 2*(floor(kymHLen/lineSegLen)-1);
speedLimit  = floor(kymHLen/lineSegLen)-1;

%Calculate a score for each discretized speed or line slope in the range
% (-speedLimit,speedLimit). The step size if one pixel.
% v     : flow speed or the slope of lines in the kymograph.
% score : Correlation score for each line slope.
score  = zeros(2*speedLimit+1,1);
%Negative speed.
for v = -speedLimit:-1
   k = v+speedLimit+1;

   numLineSegs = 0;
   %We calculate scores for all line segments of length 'lineSegLen' and in
   % the direction of slope, v.
   %First, check all the lines that start from each pixel on the top edge of 
   % the banded kymograph, 'bandKym'.
   for j = 2:kymHLen
      %The total length of the line in terms of frames.
      totalLen = min(numFrames,floor((j-1)/-v)+1);
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
         lineScore = lineSegLen* ...
            prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))))^(2/lineSegLen);
         if lineScore >= lineThreshold
            score(k) = score(k) + lineScore;
         end
      end
   end

   %Then, check all the lines that start from each pixel on the right edge of 
   % the banded kymograph, 'bandKym'.
   for j = 2:numFrames
      totalLen = min(numFrames-j+1,floor((kymHLen-1)/-v)+1);
      for f2 = lineSegLen:totalLen
         numLineSegs = numLineSegs+1;
         f1 = f2-lineSegLen+1;
         f  = f1:f2;

         col = kymHLen+(f-1)*v;
         ind = (col-1)*numFrames+f+j-1;

         lineScore = lineSegLen* ...
            prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))))^(2/lineSegLen);
         if lineScore >= lineThreshold
            score(k) = score(k) + lineScore;
         end
      end
   end

   if numLineSegs ~= 0
      score(k) = score(k)/numLineSegs;
   end
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

      lineScore = lineSegLen* ...
         prod(bandKym(f,j)/sqrt(sum(bandKymP2(f,j))))^(2/lineSegLen);
      if lineScore >= lineThreshold
         score(k) = score(k) + lineScore;
      end
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
      totalLen = min(numFrames,floor((j-1)/v)+1);
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
         lineScore = lineSegLen* ...
            prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))))^(2/lineSegLen);
         if lineScore >= lineThreshold
            score(k) = score(k) + lineScore;
         end
      end
   end

   %Then, check all the lines that start from each pixel on the left edge of 
   % the banded kymograph, 'bandKym'.
   for j = 2:numFrames
      totalLen = min(numFrames-j+1,floor((kymHLen-1)/v)+1);
      for f2 = lineSegLen:totalLen
         numLineSegs = numLineSegs+1;
         f1 = f2-lineSegLen+1;
         f  = f1:f2;

         col = 1+(f-1)*v;
         ind = (col-1)*numFrames+f+j-1;

         lineScore = lineSegLen* ...
            prod(bandKym(ind)/sqrt(sum(bandKymP2(ind))))^(2/lineSegLen);
         if lineScore >= lineThreshold
            score(k) = score(k) + lineScore;
         end
      end
   end

   if numLineSegs ~= 0
      score(k) = score(k)/numLineSegs;
   end
end

vv = [1:length(score)]-speedLimit-1;

%Get rid of zero scores. These are those speeds that are out of the range of
% detection when the speedLimit is set too high.
ind = find(score<0.1);
score(ind) = [];
vv(ind)    = [];

%We use least-square B-spline to interpolate the score samples with coarser
% breaks and find all local maximum. We then do a significance test where the
% difference between the local maximum and the mean has to be significantly 
% bigger than the deviation. All maximums that pass the significance test will
% have a consistence test. The local maximum around which the variance between
% the sampled score and the B-spline interpolation is the smallest will be the
% slope (or speed) we look for.
brk  = linspace(vv(1),vv(end),max(6,floor(length(vv)/6)));
knot = augknt(brk,4);
sp   = spap2(knot,4,vv,score.');
dsp  = fnder(sp);

%To make sure we find all the local maximum, we go through consecutive
% intervals of every 3 sampling score points.
locMax  = [];
locMaxS = [];
lmCount = 0; %local maximum count.
for k = 3:2:length(vv)-4
   [v,maxScore] = fminbnd(@scoreFun,vv(k-1),vv(k+1),[],sp);
   if v > vv(1) & v < vv(end) & abs(fnval(dsp,v)) < 1e-4
      %True local maximum
      lmCount          = lmCount+1;
      locMax(lmCount)  = v;
      locMaxS(lmCount) = -maxScore;
   end
end
%The last interval.
[v,maxScore] = fminbnd(@scoreFun,vv(k-1),vv(end-1),[],sp);
if v > vv(1) & v < vv(end) & abs(fnval(dsp,v)) < 1e-4
   %True local maximum
   lmCount          = lmCount+1;
   locMax(lmCount)  = v;
   locMaxS(lmCount) = maxScore;
end

%Compute mean and standard deviation. When compute the deviation, we only
% compute the deviation of scores higher than the mean.
scoreInterp = fnval(sp,vv(2:end-1));
avg   = sum(scoreInterp)/length(scoreInterp);
hsInd = find(scoreInterp>=avg);
dev   = sum(abs(scoreInterp(hsInd)-avg))/length(hsInd);

%Significance test.
sigThreshold = 1.2;
inSigI       = [];
for k = 1:length(locMax)
   if locMaxS(k)-avg < sigThreshold*dev
      inSigI = [inSigI k];
   end
end
locMax(inSigI)  = [];
locMaxS(inSigI) = [];

%Consistence test.
if isempty(locMax)
   %No scores are significant enough to be trusted as the speed. It means 
   % the kymograph is not clear enough to distinguished by our score function.
   v = 0;
   return;
else
   locVar = zeros(size(locMax));
   for k = 1:length(locMax)
      %Get the closest sampling point to the local maximum.
      %v         = floor(locMax(k)+0.5);
      %vI        = find(vv==v); %The corresponding index into the score.
      [tmp,ind] = min(abs(vv-v));
      vI        = [max(1,ind-2):min(length(vv),ind+2)];
      locVar(k) = sum((score(vI).'-fnval(sp,vv(vI))).^2);
   end
end

%Choose the one with minimum local variance between sampling scores and
% interpolation.
[minV,k] = min(locVar);
v        = locMax(k);
%ind = find((locMaxS-avg)>0.9*(locMaxS(k)-avg));
%[v,k] = min(abs(locMax(ind)));
%v = locMax(k);

%We remove max scores that may give false speed. A max score is trusted only
% when the differences between its nearest two neighbers and the average are 
% higher than the deviation.
%sigScoreIndex      = 1:length(score); %Significant scores.
%minNumScoreSamples = max(5,ceil(length(score)/2));
%trueMax            = 0; %First assume the max is false.
%while trueMax == 0 & length(sigScoreIndex) > minNumScoreSamples
%   [maxScore,k] = max(score(sigScoreIndex));
%
%   %if k > (length(score)+1)/2 + 2
%   %   [maxScore2,k2] = max(score((end+1)/2:k-2));
%   %   k2 = k2+(length(score)+1)/2-1;
%   %   if k2 ~= k-2
%   %      maxScore = maxScore2;
%   %      k = k2;
%   %   end
%   %elseif k < (length(score)-1)/2 - 2
%   %   [maxScore2,k2] = max(score(k+2:(end+1)/2));
%   %   k2 = k2+k+1;
%   %   if k2 ~= k+2
%   %      maxScore = maxScore2;
%   %      k = k2;
%   %   end
%   %end
%
%   %Check if the difference between the max score and the average is 
%   % significant enough to accept the identified speed.
%   avg = mean(score(sigScoreIndex));
%   hsIndex = find((score(sigScoreIndex)-avg)>=0);
%   dev = sqrt(sum((score(sigScoreIndex(hsIndex))-avg).^2)/length(hsIndex));
%   %dev = std(score(sigScoreIndex),1);
%   if k == 1
%      testI = [k k+1 k+2];
%   elseif k == length(sigScoreIndex)
%      testI = [k-2 k-1 k];
%   else
%      testI = [k-1 k k+1];
%   end
%
%   if min(score(sigScoreIndex(testI))-avg) < dev
%      %This max score failed the test and is removed.
%      sigScoreIndex(k) = [];
%   else
%      trueMax = 1;
%   end
%end
%
%if trueMax == 0
%   v = 0;
%   return;
%end
%
%v  = sigScoreIndex(k)-speedLimit-1;
%Fine tune: subpixel velocity tracking. We first interpolate the speed at
% sampling slopes with cubic spline and then find the maximum.
%brk  = linspace(vv(1),vv(end),max(5,floor(length(vv)/2)));
%knot = augknt(brk,4);
%sp   = spap2(knot,4,vv,score.');
%%pp = csape(vv,score);
%v  = fminbnd(@scoreFun,max(vv(1),v-1),min(vv(end),v+1),[],sp);

%Since we interpolate the score to half pixels, we have to divide v by 2 to
% get the real speed.
v = v/2;


function f = scoreFun(v,pp)

f = -fnval(pp,v);
