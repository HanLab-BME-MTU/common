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
%                vector (counter clockwise). 
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

numFineDir   = 6;
numCoarseDir = 40;
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
theta    = linspace(0,pi,numCoarseDir);
smplStep = theta(2)-theta(1);
thetaExt = [theta(end-4:end-1)-pi theta(1:end-1) theta(1:5)+pi];

%Set up knots sequence for least-square interpolation.
numBreaks = max(8,ceil(length(thetaExt)/4));
brk       = linspace(thetaExt(1),thetaExt(end),numBreaks);
knot      = augknt(brk,4);
brkI      = brk(2)-brk(1);

numPoints = length(x);
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
         waitbar(j/numTotalDir);
      end
   end

   %Use least-square interpolation.
   vExt     = [-v(end-4:end-1);v(1:end-1);-v(1:5)].';

   %Get rid of some bad direction where the calculated speed from kymograph is
   % not realistic. Unrealistic means that the speed is 3 times the average.
   avgV = sum(abs(vExt))/length(vExt);
   goodVI = find(abs(vExt)<3*avgV);
   goodTheta = thetaExt(goodVI);
   goodVel   = vExt(goodVI);
   sp        = spap2(knot,4,goodTheta,goodVel);
   [locThetaM,locSpeedM] = calLocMax(sp,[thetaExt(3:3:end-4) thetaExt(end-2)]);

   [maxSpeed,ind] = max(locSpeedM);
   thetaMax       = locThetaM(ind);

   %Since 'maxSpeed' is from least-square B-spline interpolation, it might be
   % biased from the sampled highest speed. We check in the interval
   % [thetaM-brkI,thetaM+brkI] to get the highest sampling speed.
   ind = find(abs(goodTheta-thetaMax)<brkI);
   [tmp,locInd] = max(abs(goodVel(ind)));
   maxInd       = goodVI(ind(locInd));
   %[maxSpeed,ind] = max(abs(v));

   %Fine tune by adding more sampling directions around the current max speed 
   % direction.
   leftI = max(1,maxInd-3);
   rightI = min(length(thetaExt),maxInd+3);
   thetaF = (thetaExt(leftI:rightI-1)+thetaExt(leftI+1:rightI))/2;
   for j = 1:length(thetaF)
      lineX = x(k)+ len*cos(thetaF(j))/2*[-1 1]; 
      lineY = y(k)+ len*sin(thetaF(j))/2*[-1 1]; 

      %Generate the kymograph along the line given by 'lineX' and 'lineY'.
      kym   = imKymograph(stack,lineX,lineY,width,'verbose','off');
      vF(j) = imKymoSpeed(kym,width);

      if strcmp(verbose,'on') == 1
         waitbar((j+numCoarseDir)/numTotalDir);
      end
   end

   %Get rid of bad directions.
   goodVFI = find(vF<3*avgV);
   goodThetaF = thetaF(goodVFI);
   goodVF     = vF(goodVFI);
   brk  = thetaExt(leftI:rightI);
   knot = augknt(brk,4);
   sp   = spap2(knot,4,[goodTheta goodThetaF],[goodVel goodVF]);
   [locThetaM,locSpeedM] = calLocMax(sp,brk);

   [maxSpeed,ind] = max(locSpeedM);
   thetaMax       = locThetaM(ind);

   %[tmp,ind] = min(abs(thetaExt-thetaMax));
   %thetaF    = [(thetaExt(ind-2)+thetaExt(ind-1))/2 ...
   %             thetaExt(ind-1)+smplStep/3 ...
   %             thetaExt(ind)-smplStep/3 thetaExt(ind)+smplStep/3 ...
   %             thetaExt(ind+1)-smplStep/3 (thetaExt(ind+2)+thetaExt(ind+1))/2];
   %vF        = zeros(size(thetaF));
   %for j = 1:numFineDir
   %   %Coordinates of the two end points of the line throught point
   %   % (x(k),y(k)).
   %   lineX = x(k)+ len*cos(thetaF(j))/2*[-1 1]; 
   %   lineY = y(k)+ len*sin(thetaF(j))/2*[-1 1]; 

   %   %Generate the kymograph along the line given by 'lineX' and 'lineY'.
   %   kym   = imKymograph(stack,lineX,lineY,width,'verbose','off');
   %   vF(j) = imKymoSpeed(kym,width);

   %   if strcmp(verbose,'on') == 1
   %      waitbar((j+numCoarseDir)/numTotalDir);
   %   end
   %end

   %Fine tune by interpolation.
   %Get rid of bad directions.
   %goodVFI = find(vF<3*avgV);
   %sp      = spap2(knot,4,[thetaExt(goodVI) thetaF(goodVFI)], ...
   %   [vExt(goodVI) vF(goodVFI)]);
   %[thetaMax,maxSpeed] = fminbnd(@vFun,thetaMax-smplStep, ...
   %   thetaMax+smplStep,[],sp);
   %maxSpeed = fnval(sp,thetaMax);
   %   ,[],pp);
   %if v(ind) < 0
   %   thetaMax = pi + thetaMax;
   %end

   %[maxSpeed,ind] = max(abs(vF));
   %if vF(ind) < 0
   %   thetaMax = pi + thetaF(ind);
   %end

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
%Calculate local maximum for each interval given in theta. To qualify for a
% local maximum, the derivative has to be zero.

dTheta = min(diff(theta))/5;
jj = 0;
for j = 1:length(theta)-1
   [thetaM,speedM] = fminbnd(@vFun,theta(j),theta(j+1),[],sp);
   if fnval(sp,thetaM-dTheta) > speedM & fnval(sp,thetaM+dTheta) > speedM
      jj = jj+1;
      locThetaM(jj) = thetaM;
      locSpeedM(jj) = -speedM;
   end
end



function f = vFun(theta,pp)

f = -abs(fnval(pp,theta));
