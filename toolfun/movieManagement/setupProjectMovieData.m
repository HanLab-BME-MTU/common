function setupProjectMovieData(projectFolder,pixSize,timeInterval,zSpacing,forceReplace)
%SETUPPROJECTMOVIEDATA sets up the movieData for a collection of movies with the same properties.
%
% setupProjectMovieData(projectFolder,pixSize,timeInterval,forceReplace)
%
% This function takes movie folders (which must be set up as
% done in setMovieFolders.m) and creates their movieData using the same
% pixel size and time interval on each.
%
%
% Input:
%   projectFolder - this is the folder which contains all the movie folders.
%   This is the same folder you specified when using setupMovieFolders.m
%
%   pixSize - image pixel size in nanometers for all movies in folder
%
%   timeInterval - time interval between images in seconds for all movies in
%   folder
%   
%   zSpacing - Spacing between z-slcices, ** if a 3D movie ** If not input,
%   movie's are assumed to be 2D.
%
%   forceReplace - If true, a new movieData will be created even if there
%   is an existing one in the same directory. This will erase all
%   processing that has been logged in the movieData.
%   Optional. Default is false.
%
% 
% Output: 
% 
%   The newly created movieData structures will be saved in each movie's
%   analysis directory as a file named movieData.mat.
% 
% 
% Hunter Elliott 
% Re-written 7/2010
%


if nargin < 1 || isempty(projectFolder)
    error('You must specify a project directory, a time interval and a pixel size!')
end

if nargin < 2
    pixSize = [];
end

if nargin < 3
    timeInterval = [];
end

if nargin < 4
    zSpacing = [];
end

if nargin < 5 || isempty(forceReplace)
    forceReplace = false;
end

%Get the folders for each movie
movieFolders = dir([projectFolder filesep]);
movieFolders = movieFolders(arrayfun(@(x)(x.isdir && ... %Retain only the directories. Do it this way so it works on linux and PC
    ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),movieFolders)); 

   
nMovies = length(movieFolders);


%Go through each and try to set up the movieData
for j = 1:nMovies
            
    disp(['Setting up movie ' num2str(j) ' of ' num2str(nMovies)])           
    
    %Switch to folder
    cd([projectFolder filesep movieFolders(j).name]);
    
    clear chans;
    
    %Check for an existing movieData
    if forceReplace || ~exist([pwd 'movieData.mat'],'file')
    
        %Set this as the output directory
        outDir = pwd;    

        %Look for sub-folder named "images"        
        if exist([pwd filesep 'images'],'dir')
            
                        
            %Check for sub directories containing the images for each
            %channel
            chanDir = dir([pwd filesep 'images']);
            chanDir = chanDir(arrayfun(@(x)(x.isdir && ... %Retain only the directories.
                            ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),chanDir));
            
            %Check if it is a single- or multiple-channel movie and set up the channels            
            if ~isempty(chanDir)
            
                %Determine the number of images and image size in each channel
                %directory.
                nChanDir = numel(chanDir);
                nImages = zeros(nChanDir,1);
                imSize = zeros(2,nChanDir);                
                
                for i = 1:nChanDir

                    %Find images in this directory
                    chans(i) = Channel([pwd filesep 'images' filesep chanDir(i).name]);
                    imFiles = imDir(chans(i).channelPath_);                    
                    nImages(i) = length(imFiles);                    
                                        
                    if nImages(i) > 0
                        
                        %Load the header for one image to determine size
                        imInfo = imfinfo([chans(i).channelPath_ filesep imFiles(1).name]);
                        imSize(1,i) = imInfo.Height;
                        imSize(2,i) = imInfo.Width;
                                                
                    end

                end                                                
                
                
            else %If it's a single-channel movie...
                
                chans = Channel([pwd filesep 'images']);
                imFiles = imDir(chans(1).channelPath_);
                nImages = length(imFiles);
                if nImages > 0
                    %Load the header for one image to determine size
                    imInfo = imfinfo([chans(1).channelPath_ filesep imFiles(1).name]);
                    imSize(1,1) = imInfo.Height;
                    imSize(2,1) = imInfo.Width;                                      
                end                                
            end
            
            %Make sure at least one image was found!
            if all(nImages>=1) && numel(unique(nImages)) == 1 && size(unique(imSize),2) == 1
                %If everything worked out, create the new movieData.
                
                if isempty(zSpacing)
                                    
                    MD = MovieData(chans,outDir,outDir,'movieData.mat',[],pixSize,timeInterval);
                else
                    MD = MovieData3D(chans,outDir,pixSize,timeInterval,zSpacing);                    
                end
                
                try
                    MD.sanityCheck
                    MD.saveMovieData;                
                    disp('MovieData setup!')                                    
                catch errMess
                    %Return to original folder
                    cd(projectFolder);
                    error(['Problem setting up movie : ' errMess.message]);                    
                end
                    
                
            else
               disp(['All image sub-directories must contain the same number of images of the same size! - cannot setup movieData for folder ' pwd '!' ])
            end            
        else
           disp(['Movie folder ' pwd ' has no sub-directory named "images" - cannot setup movieData!']) 
        end
    else
        disp(['Movie ' num2str(j) ' of ' num2str(nMovies) ' already had a movieData.mat! Doing nothing...'])
        
    end
end

%Return to original folder
cd(projectFolder);