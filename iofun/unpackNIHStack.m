function unpackNIHStack()
%UNPACKNIHSTACK unpacks a tiff stack according to the definition in NIH image
%
% SYNOPSIS unpackNIHStack
%
% INPUT none;  the filename is set interactively
%
% OUTPUT none; the frames are stored in the same directory
%
% NOTE the outcoming filenames are composed by the root of the 
%      stack filename, an account of the frame and the extension '.tif'

[fName,dirName] = uigetfile('*.tif','unpackNIHStack ...');

if( isa(fName,'char') & isa(dirName,'char'))
	   info = imfinfo([dirName,fName]);
else 
   return;
end;

nImgs = size(info,2);

if(nImgs<=1)
   warning('specified file is not a stack');
   return;
end;

% get the filename
[fPath,fBody,fNo,fExt]=getfilenamebody([dirName,fName]);

if(~isempty(fNo))
   warning('stack filenames with alphanumeric characters not supported');
   return;
end;

for i=1:nImgs
   imName= strcat(fPath,filesep,fBody,int2str(i),fExt);
   data = imread([dirName,fName],i);
   imwrite(data,imName,'Compression','none');
end;