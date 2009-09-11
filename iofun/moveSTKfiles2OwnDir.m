function moveSTKfiles2OwnDir

%allow user to choose directory
basePath = uigetdir([],'Please select directory of STK files of interest');

%go to the chosen directory
cd(basePath);

%find all files ending with ".nd" in the chosen directory
fileList = searchFiles('.nd',[],basePath,1);
numList = size(fileList,1);

%go over all movies
for iFile = numList : -1 : 1
   
    %get the name of the movie
    movieName = fileList{iFile,1};
    movieName = movieName(1:end-3);

    %generate new directory name based on movie name
    newDir = fullfile(basePath,movieName);

    %make the new directory
    mkdir(fullfile(basePath,movieName));
    
    %go to the new directory
    cd(newDir);
    
    %move all the files belonging to this movie into its directory
    movefile([newDir '*']);
    
    %go back to original directory
    cd(basePath);
    
end
