function mismatchRgbMergeBatchSetup(inputDir)
%This program (will) read through a folder containing two exposures of
%vimentin in folders with 'Dim' and 'Bright' in their title and loop
%through the folder to create hdrMerge images in their own folder for each
%one.
%-Jessica Tytell December 5, 2011


%Get the folders for each movie
movieFolders = dir([inputDir filesep '*EB*VIM*']);
movieFolders = movieFolders(arrayfun(@(x)(x.isdir && ... %Retain only the directories. Do it this way so it works on linux and PC
    ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),movieFolders));


nMovies = length(movieFolders);

%loop through movies
for j = 1:nMovies
    
    disp(['Processing folder ' num2str(j) ' of ' num2str(nMovies)])
    
    %Get current folder path for readability
    currDir = [inputDir filesep movieFolders(j).name];
    disp(currDir);
    
    %get bright and dim file directories 
    baseStruct = dir([currDir filesep 'images']);
    intermitStruct = dir([currDir filesep 'hdrMerge']);
    disp(baseStruct);
    disp(intermitStruct);
    
    %get directory names to pass to batch function
    if length(baseStruct) > 1
        disp('Too many folders labeled "images" ');
    else 
        dimName = dimStruct.name;
        dimDir = [currDir filesep dimName];
    end
    
    if length(brightStruct) > 1
        disp('Too many folders labeled "Bright"');
    else
        brightName = brightStruct.name;
        brightDir = [currDir filesep brightName];
    end
    
    %send to hdrMergeFileBatch
    hdrMergeFileBatch(dimDir, brightDir);
    
end
