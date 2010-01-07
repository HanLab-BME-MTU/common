function image = stkRead(filename)
% 
% image = stkRead(filename)
% 
% This function is just a wrapper function for metaTiffRead which
% simplifies it's use by automatically loading all images in a stack and
% concatenating them into a 3D matrix, throwing away all the tag
% information and suppressing command prompt output.
% 
% 
% The output will be a 3d matrix containing the image values.
%
% Hunter Elliott, 11/2009
%

image = metaTiffRead(filename,[],[],0);

image = cat(3,image(:).data);
