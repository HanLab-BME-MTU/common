function [vx,vy]= imKymoVelocity(stack,x,y,len,width,varargin)
%imKymoVelocity: Calculate the flow velocity from a stack of movie images. 
%
% SYNOPSIS :
%    [vx,vy] = imKymoVelocity(stack,x,y,len,width)
%    [vx,vy] = imKymoVelocity(stack,x,y,len,width,varargin)
%    [vx,vy] = imKymoVelocity(stack,x,y,len,width,'verbose','on')
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
%
% OUTPUT :
%    vx : x-component of the velocity vector.
%    vy : y-component of the velocity vector.

if mod(nargin,2) == 0
   error('The optional parameter/value has to appear in pair.');
end

numDirections = 50;

verbose = 'on';

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
vx = zeros(size(x));
vy = zeros(size(x));
numPoints = length(x);
for k = 1:numPoints
   v     = zeros(numDirections,1);
   theta = linspace(0,pi,numDirections);

   if strcmp(verbose,'on') == 1
      wbH = waitbar(0,['Processing point ' num2str(k) ' out of ' ...
         num2str(numPoints) ' ...']);
      for j = 1:numDirections
         %Coordinates of the two end points of the line throught point
         % (x(k),y(k)).
         lineX = x(k)+ len*cos(theta(j))/2*[-1 1]; 
         lineY = y(k)+ len*sin(theta(j))/2*[-1 1]; 

         %Generate the kymograph along the line given by 'lineX' and 'lineY'.
         kym{j} = imKymograph(stack,lineX,lineY,width,'verbose','off');
         v(j)   = imKymoSpeed(kym{j},width);
         waitbar(j/numDirections);
      end
   else
      for j = 1:numDirections
         %Coordinates of the two end points of the line throught point
         % (x(k),y(k)).
         lineX = x(k)+ len*cos(theta(j))/2*[-1 1]; 
         lineY = y(k)+ len*sin(theta(j))/2*[-1 1]; 

         %Generate the kymograph along the line given by 'lineX' and 'lineY'.
         kym{j} = imKymograph(stack,lineX,lineY,width,'verbose','off');
         v(j)   = imKymoSpeed(kym{j},width);
      end
   end

   [maxSpeed,ind] = max(abs(v));
   vx(k) = v(ind)*cos(theta(ind));
   vy(k) = v(ind)*sin(theta(ind));
end

if strcmp(verbose,'on') == 1
   close(wbH);
end
