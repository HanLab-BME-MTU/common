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

%The length of the line (or curve) in the movie image. It is also the
% horizontal length of the kymograph.
kymHLen = size(kym,2);

%Scale the intensity of the kymograph to be from 0.1 to 1 to increase contrast
% and to increase the difference between scores.
kym     = double(kym);
minKymI = min(kym(:));
maxKymI = max(kym(:));

if maxKymI-minKymI ~= 0
   kym = (kym-minKymI)/(maxKymI-minKymI)*1 + 0.1;
end

%We choose the maximum speed so that at this speed, a speckle will 
% flow 1/3 (if lineSeg = 3) the whole length of the selected curve for
% kymograph analysis.
speedLimit  = floor(kymHLen/2);

%Calculate the score function of speed.
[score,vv] = calScore(kym,numFrames,[-speedLimit,speedLimit],0.5);

%We use least-square B-spline to interpolate the score samples with coarser
% breaks and test the quality of the kymograph.
numBrks = max(6,floor(length(vv)/4));
if length(vv) < 2*(numBrks-1)
   sp = csape(vv,score.');
   brk  = vv;
   brkI = (vv(end)-vv(1))/(length(vv)-1);
else
   brk  = linspace(vv(1),vv(end),numBrks);
   knot = augknt(brk,4);
   sp   = spap2(knot,4,vv,score.');
   brkI = brk(2)-brk(1);
end

%Compute mean and standard deviation of the B-spline to be used for quality
% test. When we compute the deviation, we only consider scores higher than 
% the mean.
scoreInterp = fnval(sp,vv).';
avg   = sum(scoreInterp)/length(scoreInterp);
hsInd = find(scoreInterp>=avg);
dev   = sum(abs(scoreInterp(hsInd)-avg))/length(hsInd);

%When the deviation is too small (the threshold is an empirical value)
%smDevThreshold = 1e-3;
%if dev < smDevThreshold
%   v = 0;
%   return;
%end

%find all local maximum. We then do a significance test where the
% difference between the local maximum and the mean has to be significantly 
% bigger than the deviation. 
[locMax,locMaxS] = calLocMax(sp,brk);

%Significance test.
sigThreshold = 1.5;
inSigI       = [];
for k = 1:length(locMax)
   if locMaxS(k)-avg < sigThreshold*dev
      inSigI = [inSigI k];
   end
end
locMax(inSigI)  = [];
locMaxS(inSigI) = [];

if isempty(locMax) | length(locMax) > 2
   %No scores are significant enough to be trusted as the speed or too many
   % local maximum. It means the kymograph is not clear enough to 
   % distinguished by our score function.
   v = 0;
   return;
end

if length(locMax) == 2 
   if abs(locMax(2)-locMax(1)) > brkI/2
      v = 0;
      return;
   else
      [locMaxS,ind] = max(locMaxS);
      locMax        = locMax(ind);
   end
end

inSigI       = [];
for k = 1:length(locMax)
   if abs(locMax(k)) > 1
      %We also require that the average of 'score' on [locMax(k) 0] is
      % significantly greater than the average on [0 locMax(k)].
      numSmps = 20;
      smpL = linspace(locMax(k),0,numSmps);
      smpR = linspace(0,-locMax(k),numSmps);
      avgL = sum(fnval(sp,smpL))/numSmps;
      avgR = sum(fnval(sp,smpR))/numSmps;

      if avgL-avgR < dev*abs(locMax(k))/abs(vv(1));
         inSigI = [inSigI k];
      end
   end
end
locMax(inSigI)  = [];
locMaxS(inSigI) = [];

%Oscillatory test: we calculate the difference betwee the sample scores and
% the spline. If the difference is not significantly less than the difference
% between the maximum score and the mean, the quality of the kymograph is not
% good enough and zero speed is returned.
intDiff = sum(score-scoreInterp)/length(vv);
if max(locMaxS) < intDiff*sigThreshold
   v = 0;
   return;
end

%Refinement and Consistence test. In this step, we interpolate the kymograph
% to 0.25 pixels around each local maximum (2 pixels to the left and right) 
% and calculate the score function of
% sample speeds with 0.1 pixels difference. And after least square B-spline 
% interpolation, the global maximum score will be the final candidate.
for k = 1:length(locMax)
   %Since 'locMax' is from least-square B-spline interpolation, it might be
   %biased from the sampled highest score. Here we get the speed with
   % the highest local sampling score. We check in the interval
   % [localMax(k)-brkI,localMax(k)+brkI]. 
   ind    = find(abs(vv-locMax(k))<=brkI);
   [tmp,locInd] = max(score(ind)); 
   maxInd = ind(locInd);

   [locScore,locV] = calScore(kym,numFrames,[vv(maxInd)-2,vv(maxInd)+2],0.25);
   %brk  = linspace(locV(1),locV(end),ceil(length(locV)/2));
   locBrk  = [locV(1:2:end-3) locV(end)];
   locKnt = augknt(locBrk,4);
   locSp  = spap2(locKnt,4,locV,locScore.');

   %update local maximum.
   [loclocMax,loclocMaxS] = calLocMax(locSp,locBrk);
   [locMaxS(k),locInd]    = max(loclocMaxS);
   locMax(k)              = loclocMax(locInd);
end

[maxS,maxInd] = max(locMaxS);
v = locMax(maxInd);

%%%%%%%%%%%%%%% subfunction %%%%%%%%%%%%%%%%%%%%%%%%
function f = scoreFun(v,pp)
%For the purpose of 'fminbnd'.

f = -fnval(pp,v);


%%%%%%%%%%%%%%%% subfunction %%%%%%%%%%%%%%%%%%%%%%%%%
function [locMax,locMaxS] = calLocMax(sp,vv)
%Calculate local maximum for each interval given in 'vv'. To qualify for a
% local maximum, both the left and right near neighbor has to be less.
dv = min(diff(vv))/5;

lmCount = 0; %local maximum count.
for k = 1:length(vv)-1
   [v,mVal] = fminbnd(@scoreFun,vv(k),vv(k+1),[],sp);
   if abs(fnval(sp,v-dv)) < abs(mVal) & abs(fnval(sp,v+dv)) < abs(mVal)
      %True local maximum
      lmCount          = lmCount+1;
      locMax(lmCount)  = v;
      locMaxS(lmCount) = abs(mVal);
   end
end


%%%%%%%%%%%%%%%%%%% subfunction %%%%%%%%%%%%%%%%%%%%%%%
function [score,vv] = calScore(kym,numFrames,spdRng,spdStep)
%Calculate a score for each discretized speed or line slope in the speed
% range.
% INPUT :
%    kym       : The kymograph.
%    numFrames : Number of frames.
%    spdRng    : A two elements vector that specifies the low and high speeds.
%                The speed can be negative where the sign indicates the
%                direction of the flow along the kymograph line.
%                spRng(2) >= spdRng(1).
%    spdStep   : The step size of sampling speeds. 0 < spdStep <= 1.
%
% OUTPUT :
%    vv    : Sampled flow speeds or the slope of lines in the kymograph.
%    score : Correlation score for each speed or line slope.

bandWidth = size(kym,1)/numFrames;

%Interpolate to 'spdStep' pixels so that we can sample fractional pixel 
% velocity. After interpolation, we treat 'spdStep' pixels as integers. 
% In the end, the slope has to be multiplied by 'spdStep' to get the real 
% velocity.
if spdStep < 1 & spdStep > 0
   kymHLen = size(kym,2);
   ppKym   = csape([1:kymHLen],kym);
   kym = fnval(ppKym,[1:spdStep:kymHLen]);

   %'spdStep' is treated as 1 pixel after interpolation. So, we have to
   % convert speed range.
   lowSpd  = floor(spdRng(1)/spdStep);
   highSpd = ceil(spdRng(2)/spdStep);
elseif spdStep == 1
   lowSpd  = floor(spdRng(1));
   highSpd = ceil(spdRng(2));
else
   error('''spdStep'' is in the interval (0,1].');
end

if highSpd < lowSpd
   error('The speed range is not correctly set.');
end

%The square of the intensity matrix for the use of normalization later.
kymP2 = kym.*kym;

kymHLen   = size(kym,2);
score     = zeros(highSpd-lowSpd+1,1);

%Negative speed.
for v = lowSpd:min(-1,highSpd)
   k = v-lowSpd+1;

   %We calculate correlation scores every two consecutive frames in the 
   % direction of slope, v., then normalize, sum and average.
   % 
   %The maximum width of the correlating intensity bands (number of pixels in
   % the horizontal direction) that fall into the
   % scope of the kymograph depends on the slope, v.
   corrW = kymHLen+v; %'v' is negative.
   for f = 2:numFrames
      %The row index corresponding to frame 'f'.
      row = bandWidth*(f-1);

      %The correlation score between frame 'f-1' and 'f'. 
      % Note: 'v' is negative.
      corrS = sum(sum(kym(row+1:row+bandWidth,1:corrW).* ...
         kym(row-bandWidth+1:row,1-v:kymHLen))); 

      %Normalize.
      corrS = corrS/sqrt(sum(sum(kymP2(row+1:row+bandWidth,1:corrW))))/ ...
         sqrt(sum(sum(kymP2(row-bandWidth+1:row,1-v:kymHLen))));

      %Add to the total score with weight.
      score(k) = score(k) + corrS;
   end
   score(k) = score(k)/(numFrames-1);
end

%Zero and positive speed.
for v = max(0,lowSpd):highSpd
   k = v-lowSpd+1;

   corrW = kymHLen-v; %'v' is positive.
   for f = 2:numFrames
      %The row index corresponding to frame 'f'.
      row = bandWidth*(f-1);

      %The correlation score between frame 'f-1' and 'f'. 
      % Note: 'v' is negative.
      corrS = sum(sum(kym(row+1:row+bandWidth,1+v:kymHLen).* ...
         kym(row-bandWidth+1:row,1:corrW))); 

      %Normalize.
      corrS = corrS/sqrt(sum(sum(kymP2(row+1:row+bandWidth,1+v:kymHLen))))/ ...
         sqrt(sum(sum(kymP2(row-bandWidth+1:row,1:corrW))));

      %Add to the total score with weight.
      score(k) = score(k) + corrS;
   end
   score(k) = score(k)/(numFrames-1);
end

%Shift and scale to the true speeds.
vv = spdStep*([1:length(score)]+lowSpd-1);

%Get rid of zero scores. These are those speeds that are out of the range of
% detection when the speedLimit is set too high.
ind = find(score<0.1);
score(ind) = [];
vv(ind)    = [];

