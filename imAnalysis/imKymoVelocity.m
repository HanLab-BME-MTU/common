function varargout = imKymoVelocity(stack,x,y,len,width,varargin)
%imKymoVelocity: Calculate the flow velocity from a stack of movie images. 
%
% SYNOPSIS :
%    [vx,vy] = imKymoVelocity(stack,x,y,len,width)
%    [vx,vy] = imKymoVelocity(stack,x,y,len,width,varargin)
%    [vx,vy] = imKymoVelocity(stack,x,y,len,width,'verbose','on')
%    [v,theta] = imKymoVelocity(stack,x,y,len,width,'output','angle')
%
% INPUT :
%    stack : image stack or filename of the first image or a cell array
%            of the list of files to be analyzed.
%            (image stack not yet implemented)
%            if [] a gui pops up to define the first filename  
%    x     : x-coordinates of a set of points where the velocity is
%            calculated.
%    y     : y-coordinates of a set of points where the velocity is
%            calculated.
%    len   : The length of the line along which the kymograph is produced.
%    width : The width of the band around the line along which the kymograph 
%            is produced. It must be odd integer number. If the input is even,
%            it will be added by 1.
%
%    The following optional parameters can be set:
%    'verbose' : 'on' (default) or 'off'.
%                Specify if progressing message is displayed.
%    'output'  : 'vector' (default) or 'angle'.
%                Specify if the velocity is output as components of a vector
%                or as a magnitude plus the angle between x-axis and the
%                vector (going from positive x-axis to postitive y-axis). 
%
% OUTPUT :
%    vx    : x-component of the velocity vector.
%    vy    : y-component of the velocity vector.
%    v     : Magnitude of the velocity vector.
%    theta : The angle between the x-axis and the velocity vector (counter
%            clockwise).

if mod(nargin,2) == 0
   error('The optional parameter/value has to appear in pair.');
end

numFineDir   = 8;
numCoarseDir = 20;
numTotalDir  = numFineDir+numCoarseDir;

verbose = 'on';
output  = 'vector';

%Check optional parameters.
par = [];
j   = 1;
for k = 1:2:nargin-5
   par{j}   = varargin{k};
   value{j} = varargin{k+1};
   j = j+1;
end

for j = 1:length(par)
   switch par{j}
      case 'verbose'
         if ~ischar(value{j}) | (strcmp(value{j},'on') == 0 & ...
            strcmp(value{j},'off') == 0)
            error('The value for parameter ''verbose'' is not valid.');
         end
         verbose = value{j};
      case 'output'
         if ~ischar(value{j}) | (strcmp(value{j},'vector') == 0 & ...
            strcmp(value{j},'angle') == 0)
            error('The value for parameter ''output'' is not valid.');
         end
         output = value{j};
      otherwise
         error(['Parameter' par{j} 'is not recogonized.']);
   end
end

% check whether the stack is an empty matrix
if(isempty(stack))
   [fName,dirName] = uigetfile('*.tif','imKymograph ...');
   stack = [dirName,filesep,fName];
end;

% check the width
if(mod(width,2) == 0)
   width = width + 1;
end

if length(x) ~= length(y)
   error('The lengths of the x and y-coordinate vectors are not the same.');
end


%Calculate the velocity for each point. The idea is to draw a line of length
% 'len' and width 'width' through the point in diffent directions and
% calculate the speed for each direction from the corresponding kymograph
% using the function 'imKymoSpeed'. The direction that gives the maximum speed
% is flow direction or velocity vector.
varargout{1} = zeros(size(x));
varargout{2} = zeros(size(x));

%Sampling angles.
theta = linspace(0,pi,numCoarseDir);
thetaExt = [theta(end-8:end-1)-pi theta(1:end-1) theta(1:9)+pi];

%Set up knots sequence for least-square interpolation.
numBreaks = max(5,ceil(length(thetaExt)/4));
brk       = linspace(thetaExt(1),thetaExt(end),numBreaks);
knot      = augknt(brk,4);
brkI      = brk(2)-brk(1);
brkF      = linspace(thetaExt(1),thetaExt(end),2*numBreaks);
knotF     = augknt(brkF,4);

sigThreshold  = 1.2;
osciThreshold = 1;
numPoints     = length(x);
for k = 1:numPoints
   v     = zeros(numCoarseDir,1);

   if strcmp(verbose,'on') == 1
      wbH = waitbar(0,['Processing point ' num2str(k) ' out of ' ...
         num2str(numPoints) ' ...']);
   end
   for j = 1:numCoarseDir
      %Coordinates of the two end points of the line throught point
      % (x(k),y(k)).
      lineX = x(k)+ len*cos(theta(j))/2*[-1 1]; 
      lineY = y(k)+ len*sin(theta(j))/2*[-1 1]; 

      %Generate the kymograph along the line given by 'lineX' and 'lineY'.
      kym  = imKymograph(stack,lineX,lineY,width,'verbose','off');
      v(j) = imKymoSpeed(kym,width);

      if strcmp(verbose,'on') == 1
         waitbar(j/numTotalDir,wbH);
      end
   end

   %Use least-square interpolation.
   vExt     = [-v(end-8:end-1);v(1:end-1);-v(1:9)].';

   %Get rid of some bad direction where the calculated speed from kymograph is
   % not realistic. Unrealistic means that the speed is 3 times the average.
   avgV = sum(abs(vExt))/length(vExt);
   goodVI = find(abs(vExt)<3*avgV);
   badVI  = find(abs(vExt)>=3*avgV);
   goodTheta = thetaExt(goodVI);
   goodVel   = vExt(goodVI);

   if ~isempty(badVI)
      %We don't really get rid of them but interpolate to avoid undersampling
      % when using 'spap2'.
      pp = csape(goodTheta,goodVel);
      vExt(badVI) = fnval(pp,thetaExt(badVI));
   end
   sp = spap2(knot,4,thetaExt,vExt);

   %Compute mean and standard deviation of the B-spline to be used for quality
   % test. When we compute the deviation, we only consider scores higher than 
   % the mean.
   vInterp = abs(fnval(sp,thetaExt));
   avg     = sum(vInterp)/length(vInterp);
   hvInd   = find(vInterp>=avg);
   dev     = sum(abs(vInterp(hvInd)-avg))/length(hvInd);

   [locThetaM,locSpeedM] = calLocMax(sp,[thetaExt(7:2:end-8) thetaExt(end-6)]);
   %Significance test.
   inSigI       = [];
   for j = 1:length(locThetaM)
      if locSpeedM(j)-avg < sigThreshold*dev
         inSigI = [inSigI j];
      end
   end
   locThetaM(inSigI)  = [];
   locSpeedM(inSigI) = [];

   %Oscillatory test: we calculate the difference between the sample speeds
   % and the spline. If the difference is not significantly less than the 
   % difference between the maximum speed and the mean, the quality of the 
   % kymograph is not good enough and zero speed is returned.
   inSigI = [];
   for j = 1:length(locThetaM)
      %First, find the sampling point that is closest to locThetaM(j).
      [tmp,ind] = min(abs(thetaExt-locThetaM(j)));

      %Then look to the right and left of 'ind' for indices where speeds
      % are local minimum or are less than the average speed. The sampling
      % points between the two (left and right) indices are used for
      % oscillator test.
      %To right:
      jj = ind+1;
      while jj < length(vInterp) & vInterp(jj) > avg & ...
         (vInterp(jj) > vInterp(jj-1) | vInterp(jj) > vInterp(jj+1))
         jj = jj+1;
      end
      rightI = jj-1;
      %To Left:
      jj = ind-1;
      while jj > 1 & vInterp(jj) > avg & ...
         (vInterp(jj) > vInterp(jj-1) | vInterp(jj) > vInterp(jj+1))
         jj = jj-1;
      end
      leftI   = jj+1;

      %Difference between samples and interpolation.
      intDiff = sum(abs(abs(vExt(leftI:rightI))-vInterp(leftI:rightI)))/ ...
         (rightI-leftI+1);
      if intDiff > osciThreshold*dev;
         inSigI = [inSigI j];
      end
   end
   if ~isempty(inSigI)
      locThetaM(inSigI)  = [];
      locSpeedM(inSigI) = [];
   end

   maxSpeed = 0;
   thetaMax = 0;
   while maxSpeed == 0 & ~isempty(locThetaM)
      [maxSpeed,maxI] = max(locSpeedM);
      thetaMax        = locThetaM(maxI);

      %Since 'maxSpeed' is from least-square B-spline interpolation, it might be
      % biased from the sampled highest speed. We check in the interval
      % [thetaM-brkI,thetaM+brkI] to get the highest sampling speed.
      ind = find(abs(goodTheta-thetaMax)<brkI);
      [tmp,locInd] = max(abs(goodVel(ind)));
      maxInd       = goodVI(ind(locInd));
      %[maxSpeed,ind] = max(abs(v));

      %Fine tune by adding more sampling directions around the current max
      % speed direction.
      leftI  = max(1,maxInd-2);
      rightI = min(length(thetaExt),maxInd+2);
      thetaF1 = (thetaExt(leftI:rightI-1)+2*thetaExt(leftI+1:rightI))/3;
      thetaF2 = (2*thetaExt(leftI:rightI-1)+thetaExt(leftI+1:rightI))/3;
      thetaF  = zeros(1,length(thetaF1)+length(thetaF2));
      thetaF(1:2:end) = thetaF1;
      thetaF(2:2:end) = thetaF2;
      vF     = zeros(size(thetaF));
      for j = 1:length(thetaF)
         lineX = x(k)+ len*cos(thetaF(j))/2*[-1 1]; 
         lineY = y(k)+ len*sin(thetaF(j))/2*[-1 1]; 

         %Generate the kymograph along the line given by 'lineX' and 'lineY'.
         kym   = imKymograph(stack,lineX,lineY,width,'verbose','off');
         vF(j) = imKymoSpeed(kym,width);

         if strcmp(verbose,'on') == 1
            waitbar((j+numCoarseDir)/numTotalDir,wbH);
         end
      end

      %Oscillatory test.
      intVF   = fnval(sp,thetaF);
      intDiff = sum(abs(vF-intVF))/length(thetaF);
      if intDiff <= osciThreshold*dev;
         %Get rid of bad directions.
         goodVFI    = find(vF<3*avgV);
         goodThetaF = thetaF(goodVFI);
         goodVF     = vF(goodVFI);
         spF        = spap2(knotF,4,[thetaExt goodThetaF],[vExt goodVF]);
         [locThetaMF,locSpeedMF] = calLocMax(spF,thetaExt(leftI-1:rightI+1));

         [maxSpeedF,ind] = max(locSpeedMF);
         if abs(maxSpeedF-maxSpeed) <= osciThreshold*dev
            thetaMax  = locThetaMF(ind);
            maxSpeed  = fnval(spF,thetaMax);
         else
            maxSpeed = fnval(sp,thetaMax);
         end
      else
         maxSpeed = 0;
         thetaMax = 0;
         locThetaM(maxI) = [];
         locSpeedM(maxI) = [];
      end
   end

   if strcmp(verbose,'on') == 1
      close(wbH);
   end

   if strcmp(output,'vector') == 1
      varargout{1}(k) = maxSpeed*cos(thetaMax);
      varargout{2}(k) = maxSpeed*sin(thetaMax);
   elseif strcmp(output,'angle') == 1
      varargout{1}(k) = maxSpeed;
      varargout{2}(k) = thetaMax;
   end
end

%%%%%%%%%%%% subfunction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [locThetaM,locSpeedM] = calLocMax(sp,theta)
%Calculate local maximum for each interval given in 'theta'. To qualify for a
% local maximum, both the left and right near neighbor has to be less.

dTheta = min(diff(theta))/10;

thetaM = zeros(1,length(theta)-1);
speedM = zeros(1,length(theta)-1);
lmCount = 0; %Local maximum count.
for j = 1:length(theta)-1
   [thetaM(j),speedM(j)] = fminbnd(@vFun,theta(j),theta(j+1),[],sp);
   if abs(fnval(sp,thetaM(j)-dTheta)) < abs(speedM(j)) & ...
      abs(fnval(sp,thetaM(j)+dTheta)) < abs(speedM(j))
      if lmCount == 0 | ...
         (lmCount > 0 & abs(thetaM(j)-locThetaM(lmCount)) > dTheta)
         lmCount            = lmCount+1;
         locThetaM(lmCount) = thetaM(j);
         locSpeedM(lmCount) = -speedM(j);
      else
         if abs(speedM(j)) > abs(locSpeedM(lmCount))
            locThetaM(lmCount) = thetaM(j);
            locSpeedM(lmCount) = -speedM(j);
         end
      end
   end
end

if lmCount == 0
   [locSpeedM,ind] = min(speedM);
   locThetaM       = thetaM(ind);
   locSpeedM       = -locSpeedM;
end



function f = vFun(theta,pp)

f = -abs(fnval(pp,theta));
