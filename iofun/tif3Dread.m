function image = tif3Dread(filename)

% 
% image = tif3Dread(filename)
% 
% This reads every image in the specified 3D (multi-page) .tif file and
% combines all the images into a single 3D matrix.
% 
% Input:
% 
%   filename - The name of the multi-page .tif file to read.
% 
% 
% Output:
% 
%   image - The 3D matrix containing the concatenated images.
%
%
% Hunter Elliott
% 1/2010
%

if nargin < 1 || isempty(filename)
    error('Must input a file name!')
end

if ~exist(filename,'file')
    error('Specified file does not exist!');
end

%get the indices of each page in the selected .tif file
i = 1:length(imfinfo(filename));

%Load them all
image = arrayfun(@(x)(imread(filename,x)),i,'UniformOutput',false);

%Concatenate them into a 3D array
image = cat(3,image{:});