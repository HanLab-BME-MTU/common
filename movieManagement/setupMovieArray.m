function movieArray = setupMovieArray(parentDirectory)

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

%Allow the user to select among the files found
[iSel,selectedFiles] = listSelectGUI(fList,[],'move');

%Load all of the movie data files and put them in an array
nFiles = length(selectedFiles);

movieArray = cell(1,nFiles);

for j = 1:nFiles
    
    tmp = load(selectedFiles{j},'movieData');
    movieArray{j} = tmp.movieData;
    
end
