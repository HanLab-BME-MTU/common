function setupProjectMovieData(projectFolder,pixSize,timeInterval,forceReplace)
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
% Re-written 1/20/2010
%


if nargin < 3 || isempty(pixSize) || isempty(timeInterval)
    error('You must specify a project directory, a time interval and a pixel size!')
end

if nargin < 4 || isempty(forceReplace)
    forceReplace = false;
end

%Get the folders for each movie
movieFolders = dir([projectFolder filesep]);
movieFolders = movieFolders(arrayfun(@(x)(x.isdir && ... %Retain only the directories. Do it this way so it works on linux and PC
    ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),movieFolders)); 

   
nMovies = length(movieFolders);


%Set common properties
movieData.pixelSize_nm = pixSize;
movieData.timeInterval_s = timeInterval;

%Go through each and try to set up the movieData
for j = 1:nMovies
            
    disp(['Setting up movie ' num2str(j) ' of ' num2str(nMovies)])        
    
    %Re-set fields from previous movieData
    movieData.nImages = [];
    movieData.imSize =  [];
    movieData.channelDirectory = {};
    
    
    %Switch to folder
    cd([projectFolder filesep movieFolders(j).name]);
    
    %Check for an existing movieData
    if forceReplace || ~exist([pwd 'movieData.mat'],'file')
    
        %Set this as the analysis folder
        movieData.analysisDirectory = pwd;    

        %Look for sub-folder named "images"        
        if exist([pwd filesep 'images'],'dir')
            
            %Set this as the movie's image directory
            movieData.imageDirectory = [movieData.analysisDirectory filesep 'images'];
                        
            %Check for sub directories containing the images for each
            %channel
            chanDir = dir([pwd filesep 'images']);
            chanDir = chanDir(arrayfun(@(x)(x.isdir && ... %Retain only the directories.
                            ~(strcmp(x.name,'.') || strcmp(x.name,'..'))),chanDir));
            
            %Check if it is a single- or multiple-channel movie and set up the channels            
            if ~isempty(chanDir) > 0
            
                %Determine the number of images and image size in each channel
                %directory.
                for i = 1:length(chanDir)

                    %Find images in this directory
                    imFiles = imDir([pwd filesep 'images' filesep chanDir(i).name]);

                    movieData.nImages(i) = length(imFiles);
                    movieData.channelDirectory{i} = chanDir(i).name;
                    if movieData.nImages(i) > 0
                        %Load the header for one image to determine size
                        imInfo = imfinfo([pwd filesep 'images' filesep chanDir(i).name filesep imFiles(1).name]);
                        movieData.imSize(1,i) = imInfo.Height;
                        movieData.imSize(2,i) = imInfo.Width;                
                    end


                end
            else %If it's a single-channel movie
                movieData.channelDirectory{1} = '';                
                imFiles = imDir([pwd filesep 'images']);
                movieData.nImages = length(imFiles);
                if movieData.nImages > 0
                    %Load the header for one image to determine size
                    imInfo = imfinfo([pwd filesep 'images' filesep imFiles(1).name]);
                    movieData.imSize(1,1) = imInfo.Height;
                    movieData.imSize(2,1) = imInfo.Width;                
                end                                
            end
            %Make sure at least one image was found!
            if sum(movieData.nImages(:)) > 0
                %If everything worked out, save the newly created
                %movieData.
                updateMovieData(movieData);
            else
               disp(['No images found in the image directory for movie ' pwd ' - cannot setup movieData!'])
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