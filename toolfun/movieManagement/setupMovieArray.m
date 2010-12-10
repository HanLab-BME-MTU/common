function movieArray = setupMovieArray(parentDirectory)
%SETUPMOVIEARRAY creates an array of MovieData objects by searching through a directory
% 
% movieArray = setupMovieArray;
% 
% movieArray = setupMovieArray(parentDirectory);
% 
% This function finds the MovieData object for every movie in the
% specified parent directory. These MovieData objects should be saved in
% files named MovieData.mat, and should be formatted as created by
% setupMovieDataGUI.m
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

%movieArray = MovieData(1,nFiles); Not sure how to pre-allocate without
%running into access problems to private fields...

isGood = true(nFiles,1);
for j = 1:nFiles
    
    tmp = load(selectedFiles{j});
    fNames = fieldnames(tmp);
    if numel(fNames) > 1 || ~isa(tmp.(fNames{1}),'MovieData');
        disp(['Invalid MovieData found at ' selectedFiles{j} ' - Not including in array! Please check this movieData.mat file!']);
        isGood(j) = false;
    else
        movieArray(j) = tmp.(fNames{1});
    end
    
end

movieArray = movieArray(isGood);
