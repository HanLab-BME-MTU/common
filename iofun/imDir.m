function [fileNames formatNum] = imDir(imDirectory,returnAll)
%IMDIR is a wrapper for the dir command which searches only for common image file types
% 
% fileNames = imDir(directory);
% 
% fileNames = imDir(directory,returnAll);
%
% [fileNames formatNum] = imDir(...);
%
% This function will find all files in the specified directory with common
% file extensions for images. They are returned in the same format as the
% dir command.
% 
% Input:
% 
%   directory - the directory to search for files in. (non-recursive);
%   Optional. If not input, current directory is used.
%
%   returnAll - If true, all images of any file extension will be returned.
%   If false, only the first matching set of files found will be returned,
%   in this order:
% 
%       1 - .tif 
%       2 - .TIF 
%       3 - .STK 
%       4 - .bmp 
%       5 - .BMP 
%       6 - .jpg
%       7 - .JPG
%
%   This input is optional. Default is false.
%
% Output:
%
%   fileNames - a structure containing the names and statistics for all the
%   images found. This is the same as the output of the dir function.
%
%   formatNum - the number of specified extensions in the search directory 
%
% Hunter Elliott
% 2/2010
%

%The list of supported file extensions. Feel free to add! (just update the
%help also!)
fExt = {'tif','TIF','STK','bmp','BMP','jpg'};


if nargin < 1 || isempty(imDirectory)
    imDirectory = pwd;
end

if nargin < 2 || isempty(returnAll)
    returnAll = false;
end

fileNames = [];
formatNum = 0;

for i = 1:length(fExt)
    
    tempfileNames = dir([imDirectory filesep '*.' fExt{i}]);
    if ~isempty(tempfileNames)
        formatNum = formatNum +1;
    end
    
    fileNames = vertcat(fileNames, tempfileNames);
    if ~returnAll && ~isempty(fileNames);
        break
    end
end
   


