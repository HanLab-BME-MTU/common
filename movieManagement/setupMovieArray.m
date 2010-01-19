function movieArray = setupMovieArray(parentDirectory)

% 
% movieArray = setupMovieArray;
% 
% movieArray = setupMovieArray(parentDirectory);
% 
% This function finds the movieData structure for every movie in the
% specified parent directory. These movieData structures should be saved in
% files named movieData.mat, and should be formatted as created by
% setupMovieData.m
% 
% Input:
% 
%   parentDirectory - The directory containing all the movies (and
%   movieData.mat files).
% 
% Output:
% 
%   movieArray - A cell-array containing all the movieData structures for
%   the movies in the parent directory.
%
%
%
% Hunter Elliott
% Sometime in 2008?


%Get the parent directory if not input
if nargin < 1 || isempty(parentDirectory)
    parentDirectory = uigetdir('','Select the parent directory containing all the movies:');
end

%Search for movie data files in this directory

if parentDirectory == 0 %if user clicked cancel
    movieArray = [];
    return
else
    fList = searchFiles('movieData.mat',[],parentDirectory,1,'new',1);
end

if ~isempty(fList)
    %Allow the user to select among the files found
    [iSel,selectedFiles] = listSelectGUI(fList,[],'move');
else
    error('No movieData.mat files found in specified directory! Check directory!')    
end
%Load all of the movie data files and put them in an array
nFiles = length(selectedFiles);

movieArray = cell(1,nFiles);

for j = 1:nFiles
    
    tmp = load(selectedFiles{j},'movieData');
    movieArray{j} = tmp.movieData;
    
end
