function [kym,xBand,yBand] = imKymograph(stack,x,y,width);
%IMKYMOGRAPH generates a kymograph from an image stack
% 
% SYNOPSIS kym = imKymograph(stack,x,y,width);
% 
% INPUT    stack : image stack or filename of the first image
%                  (image stack not yet implemented)
%                  if [] a gui pops up to define the first filename  
%          x     : x-coordinates of the kymographed trajectory
%          y     : y-coordinates of the kymographed trajectory
%                  The tracjectory is supposed to have no loops
%                  if length(x) = length(y) = 2 then a traditional line 
%                  kymograph is produced
%          width : width of the graph (must be odd integer number)
%
% OUPUT    kym   : image with kymograph with the dimensions 
%                  n * width x stretch(x,y)
%                  where n is the number of frames in the stack
%                  and   stretch(x,y) is the pixel length of the tracjectory

% STARTED July-9-1999 GD

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

% check whether the stack is an empty matrix
if(isempty(stack))
   [fName,dirName] = uigetfile('*.tif','imKymograph ...');
   stack = [dirName,filesep,fName];
end;

% check the width
if(mod(width,2) == 0)
   width = width + 1;
end

% check that there are no loops in the trajectory
if(length(x) ~= length(y))
   error('invalid trajectory coordinates');
end;
dx = diff(x);
dy = diff(y);
xFin(1) = x(1);
xFin(2) = x(2);
yFin(1) = y(1);
yFin(2) = y(2);
refDir = [ dx(1), dy(1)];
refDir = refDir / norm(refDir);
i = 2;
newDir = zeros(1,2);
while(i <= length(dx))
   newDir = newDir + [dx(i), dy(i)];
   newDirNorm = newDir / norm(newDir);
   if(newDirNorm * refDir' > 0)
      % the trajectory continues in forward direction
      xFin = cat(2,xFin,x(i+1));
      yFin = cat(2,yFin,y(i+1));
      refDir = newDirNorm;
      newDir = zeros(1,2);
   else
      % the trajectory continues in backward direction
   end;
   i = i+1;
end;

% get trajectory with single pixel spacing
dx = diff(xFin);
dy = diff(yFin);
s = sqrt(dx.^2+dy.^2);
sCum = cumsum(s);
sCum = [0, sCum];
sInt = 0:1:sCum(end);
xInt = interp1(sCum,xFin,sInt);
yInt = interp1(sCum,yFin,sInt);

% get coordinates of adjacent, parallel bands
pos = -(width-1)/2 : 1 : (width-1)/2;
dx = conv(xInt,[1,0,-1])/2;  % positional derivative in x
dx = dx(3:end-2);
dy = conv(yInt,[1,0,-1])/2;  % positional derivative in y
dy = dy(3:end-2);
dn = sqrt(dx.^2+dy.^2);
xBand = ones(width,1)*xInt(2:end-1) + pos' * (dy ./ dn);
yBand = ones(width,1)*yInt(2:end-1) - pos' * (dx ./ dn);

kym = [];

if(ischar(stack))
   % read image by image and do the interpolation along the band
   [fpath,fname,fno,fext]=getfilenamebody(stack);
   if(isempty(fname) | isempty(fno) | isempty(fext) )
      error('invalid first filename specified');
   end;
   
   if(~isempty(fpath))
      % change to stack directory
      oldDir = cd(fpath);
   end;
   
   dirListing = dir;
   % get all relevant filenames
   iEntry = 1;
   fileList = {};
   for( i = 1:length(dirListing))
      if(~dirListing(i).isdir)
         fileList(iEntry) = {dirListing(i).name};
         iEntry = iEntry + 1;
      end;
   end;
   
   nEntries = 0;
   imIndx = str2num(fno);
   searchName= [fname,int2str(imIndx),fext];
   while(~isempty(strmatch(searchName,fileList)))
      auxI = imread(searchName);
      % convert eventual color images to gray
      if(size(auxI,3) == 3)
         auxI = rgb2gray(auxI);
      end;
      
      % interpolate the intensities along the band 
      bandI = [];
      for(i = 1:width)
         % bandI = cat(2,bandI,improfile(auxI,xBand(i,:),yBand(i,:)));
         bandI = cat(1,bandI,interp2(double(auxI),xBand(i,:),yBand(i,:)));
      end;
      kym = cat(1,kym,bandI);
      disp(['completed file: ',searchName]);

      imIndx = imIndx + 1;
      searchName= [fname,int2str(imIndx),fext];
      
   end;
   
   % change back to original directory
   if(~isempty(oldDir))
      cd(oldDir);
   end;
else
   error('input via image stack not yet implemented');
end;

   
   


   