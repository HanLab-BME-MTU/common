function [vx,vy]= imKymoVelocity(stack,x,y,len,width)
%imKymoVelocity: Calculate the flow velocity from a stack of movie images. 
%
% SYNOPSIS :
%    [vx,vy] = imKymoVelocity(stack,x,y,len,width)
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
% OUTPUT :
%    vx : x-component of the velocity vector.
%    vy : y-component of the velocity vector.

numDirections = 50;

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
   for j = 1:numDirections
      %Coordinates of the two end points of the line throught point
      % (x(k),y(k)).
      lineX = x(k)+ len*cos(theta(j))/2*[-1 1]; 
      lineY = y(k)+ len*sin(theta(j))/2*[-1 1]; 

      %Generate the kymograph along the line given by 'lineX' and 'lineY'.
      kym  = imKymograph(stack,lineX,lineY,width,'verbose','off');
      v(j) = imKymoSpeed(kym,width);
   end
   [maxSpeed,ind] = max(abs(v));
   vx(k) = v(ind)*cos(theta(ind));
   vy(k) = v(ind)*sin(theta(ind));
end
