function movieData = setupMovieData(movieData,procedure)
%
% movieData = setupMovieData(movieData,procedure)
% 
% This movie is used to create and/or validate a movieData structure. The
% movieData structure stores data about the movie images and their
% directories, any analysis and processing that is performed on the movie,
% the parameters used in the analysis etc. 
%
%
% Input:
%
%   movieData - Optional. Structure as described above. If input, the
%   structure is checked, if not input, a new one is created and saved.
% 
%   procedure - A string describing a processing/analysis procedure to be
%   applied to the movie. A directory is created within the movie's
%   analysis directory that has the same name as the string. Also, a field
%   is created within the movieData structure to track this procedure.
% 
% 
% Output:
% 
%   movieData - The verified or newly created movieData structure. If new,
%   this has also been saved to the specified movieData's analysis
%   directory.
%
% Hunter Elliott, 2008
%

if nargin < 2
   procedure = []; 
elseif ~ischar(procedure)
    errordlg('Procedure must be a recognized string!',mfilename)
    movieData = [];
    return
end

if nargin == 0
    movieData = [];
end
   
% REMOVED THIS FEATURE AFTER THE MOVE TO BOSTON - HLE
% %If the movieData was input, covert it to the current OS. This ensures
% %linux / windows compatibility
% if ~isempty(movieData)
%     movieData = convertMovieDataOS(movieData);
% end


%Check if the analysis directory was supplied, and if not ask for it
if isempty(movieData) || ~isfield(movieData,'analysisDirectory') || ...
        isempty(movieData.analysisDirectory)   
    
    movieData.analysisDirectory = uigetdir(pwd,...
        'Select the folder for storing the analysis:');
    
    if movieData.analysisDirectory == 0
        return
    end
    
    %Check if a movieData already exists in this folder
    if exist([movieData.analysisDirectory filesep 'movieData.mat'],'file')
        dialogAns = questdlg('Previous analysis file in this folder!',...
            'Movie Analysis Setup','Use Existing','Overwrite Existing',...
            'Cancel','Cancel');
        switch dialogAns
            
            case 'Use Existing'
                tmp = load([movieData.analysisDirectory filesep 'movieData.mat'],'-mat');
                
                if ~isfield(tmp,'movieData') || ~isfield(tmp.movieData,...
                        'analysisDirectory') || isempty(tmp.movieData.analysisDirectory)
                    errordlg('Problem with existing movieData file!',mfilename)
                    movieData = [];
                    return
                else
                    movieData = tmp.movieData;
                end                
            case 'Cancel'
                movieData = [];
                return
                
        end        
    end
    
elseif ~isstruct(movieData)
    errordlg('Invalid movie data!',mfilename);
    return        
end

%Check if the image directory was supplied, and if not ask for it
if ~isfield(movieData,'imageDirectory') || ...
        isempty(movieData.imageDirectory)   
    
    movieData.imageDirectory = uigetdir(pwd,...
        'Select the folder which contains the images:');            
    
    if movieData.analysisDirectory == 0
        return
    end
    
    %Check if the folder contains images, or other folders. If it contains
    %both sub-folders and images, ask if each is a seperate channel
    folderList = dir(movieData.imageDirectory);
    imageList = imDir(movieData.imageDirectory);
    nImages = length(imageList);
    nFolders = length(folderList);
    isDir = false(nFolders,1);
    for j = 1:nFolders
        if folderList(j).isdir && ~strcmp(folderList(j).name,'.') && ~strcmp(folderList(j).name,'..')    
            isDir(j) = true;
        end
    end        
    iDir = find(isDir);
    nChannels = length(iDir);
    if nChannels > 0 
        if nImages > 0
           %Ask the user if these are channels.
           bPressed = questdlg('The image directory contains images AND sub-directories.','','Single-Channel Movie','Each sub-dir is a channel','Cancel','Cancel');
           
            switch bPressed
                
                case 'Single-Channel Movie'                    
                    movieData.nImages = nImages;                                                           
                    movieData.channelDirectory{1} = '';
                case 'Each sub-dir is a channel'
                    for j = 1:nChannels
                        movieData.channelDirectory{j} = folderList(2+j).name;
                    end 
                case 'Cancel'
                    return
                otherwise
                    return
            end                    
        else %If no images found, just assume it is multi-channel
            for j = 1:nChannels
                movieData.channelDirectory{j} = folderList(2+j).name;
            end             
        end
    else
        movieData.nImages = nImages;                                                           
        movieData.channelDirectory{1} = '';        
    end
end    
    

%Check if there is a folder for the current procedure
if ~isempty(procedure)
    
    %Setup the field for this procedure in the movieData structure
    if ~isfield(movieData,procedure)
        movieData.(procedure).status = 0;
    end
    
    %and the directory, if necessary
    if ~exist([movieData.analysisDirectory filesep procedure],'dir')
        %msgbox(['No "' procedure  '" directory - creating.'],'Movie Analysis Setup')
        mkdir([movieData.analysisDirectory filesep procedure]);
    end
    
    if ~isfield(movieData.(procedure)(1),'directory')
        movieData.(procedure)(1).directory = [ movieData.analysisDirectory filesep procedure ];
    end
end


%Check if "essential" movie information is present


%Check pixel size and get if necessary
if ~isfield(movieData,'pixelSize_nm') || isempty(movieData.pixelSize_nm)
    dialogAns = inputdlg('Enter pixel size (nm):',...
        'Pixel Size Needed',1,{'0'});    
    movieData.pixelSize_nm=str2double(dialogAns{1});
end
while (movieData.pixelSize_nm <= 0) || ~isnumeric(movieData.pixelSize_nm)...
        || ~isreal(movieData.pixelSize_nm) || isnan(movieData.pixelSize_nm)
    dialogAns = inputdlg('Enter pixel size (nm):',...
        'Pixel Size Needed',1,{'0'});
    movieData.pixelSize_nm=str2double(dialogAns{1});
end

%Check time interval and get if necessary
if ~isfield(movieData,'timeInterval_s') || isempty(movieData.timeInterval_s)
    dialogAns = inputdlg('Enter time interval (s):',...
        'Time Interval Needed',1,{'0'});    
    movieData.timeInterval_s=str2double(dialogAns{1});
end
while (movieData.timeInterval_s <= 0) || ~isnumeric(movieData.timeInterval_s)...
        || ~isreal(movieData.timeInterval_s) || isnan(movieData.timeInterval_s)
    dialogAns = inputdlg('Enter time interval (s):',...
        'Time Interval Needed',1,{'0'});
    movieData.timeInterval_s=str2double(dialogAns{1});
end

%Check number of images and get if necessary
nChan = length(movieData.channelDirectory);
if ~isfield(movieData,'nImages') || isempty(movieData.nImages)
    
    %Go through each channel and check the number of images    
    movieData.nImages = zeros(1,nChan);    
    for j = 1:nChan
        movieData.nImages(j) = length(imDir([movieData.imageDirectory filesep movieData.channelDirectory{j}]));                       
    end

end

%Check image size and get if necessary
if ~isfield(movieData,'imSize')
    
    %Go through each channel and check the size of ONE image (assume
    %they're all the same...)
    movieData.imSize = zeros(2,nChan);
    for j = 1:nChan
        currFnames = imDir([movieData.imageDirectory filesep movieData.channelDirectory{j}]);
        currIm = imread([movieData.imageDirectory filesep movieData.channelDirectory{j} filesep currFnames(j).name]);
        movieData.imSize(:,j) = size(currIm);        
    end
end


%Save the updated movie data if it is newly created. Otherwise,
%updateMovieData.m must be used to save an updated movieData
if ~exist([movieData.analysisDirectory filesep 'movieData.mat'],'file')
    save([movieData.analysisDirectory filesep 'movieData.mat'],'movieData');
end

