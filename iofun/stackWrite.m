function stackWrite(im,fileName,compType)
%STACKWRITE writes a 3D image matrix to a single multi-page .tif file 
% 
% stackWrite(im,fileName)
% stackWrite(im,fileName,compression)
% 
% This function writes the input 3D matrix to a SINGLE, multi-page .tif
% file with the specified file name. Compression is disabled by default to
% increase compatability. To write a 3D image to multiple .tif images, use
% imwritestack.m or a similar function.
% 
% Input:
%   
%   im - The 3D image matrix to write to file. The resulting image
%   bit-depth will depend on the class of this image.
% 
%   fileName - The file name to write the image to, WITH file extension.
%
%   compression - A string specifying the type of compression to use. The
%   possible options are described in the help for imwrite.m Note that many
%   of the compression types are incompatible with the stackRead.m function
%   so will require you to use tif3Dread.m to open the resulting stack.
%   Optional. Default is 'none'.
% 
% Output:
%
%   The input image matrix will be written to disk. If a file already
%   exists with that name, it will be overwritten.
%
% See Also:
%
%   stackRead.m, tif3Dread.m, imwrite.m
%
% Hunter Elliott
% 2/2010
%

if nargin < 2 || isempty(im) || isempty(fileName)
    error('You must input an image and file name!')
end

if ndims(im) ~= 3
    error('Input image must be 3D!')
end

if nargin < 3 || isempty(compType)
    compType = 'none';
end

%Write the first slice in overwrite mode, in case the file already exists.
imwrite(im(:,:,1),fileName,'tif','Compression',compType)

%All successive z-slices are appended.
for i = 2:size(im,3)        
    imwrite(im(:,:,i),fileName,'tif','WriteMode','append','Compression',compType)    
end