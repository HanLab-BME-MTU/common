function setupMovieImageFolders(projectFolder,channelNames)
%SETUPMOVIEIMAGEFOLDERS places images from different channels in their own directory based on their name
% 
% setupMovieImageFolders(projectFolder,channelNames) 
% 
% Input:
% 
% projectFolder - A directory containing all the movies to setup folders
% for. Each movie should be in it's own directory.
% 
% channelNames - The string to search for in each image name, and also the
% name of the folder it will be moved to.
% Optional. Default is basic FRET channels: CFP, FRET, DIC 
% 
% 
%Hunter Elliott 2/09
%

fExt = 'tif';%File extension of images to look for.

if nargin < 2 || isempty(channelNames)
    channelNames = {'CFP','FRET','DIC'};%default is basic FRET channels
end
nChan = length(channelNames);

%Get the folders for each movie
movieFolders = dir(projectFolder);
movieFolders = movieFolders(arrayfun(@(x)(x.isdir && ... %Retain only the directories. Do it this way so it works on linux and PC
    ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),movieFolders)); 
nMovies = length(movieFolders);


ogDir = pwd;%Store current location so we can switch back

%Go through each folder and set up the image folders
for j = 1:nMovies
        
    
    disp(['Setting up movie ' num2str(j) ' of ' num2str(nMovies)])
    
    
        
    %Switch to movie directory to save typing
    cd([projectFolder filesep movieFolders(j).name]);    
    
    %Create the images directory
    mkdir('images')
    
    
    for k = 1:nChan
    
        %Look for files fitting the current string
        imFiles = dir(['*' channelNames{k} '*.' fExt]);
        
        if ~isempty(imFiles)
            
            %Make this channel's directory
            mkdir(['images' filesep channelNames{k}]);
            
            %Move all the images into it.
            arrayfun(@(x)(movefile(x.name,['images' filesep channelNames{k} filesep x.name])),imFiles);
            
        else
            
            disp(['Couldnt find any images matching channel "' channelNames{k} '"!'])
        end
   
   
    end
    
    
end

%Switch back to the original directory.
cd(ogDir);