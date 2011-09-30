function croppedMovieData = cropMovie(movieData,cropROI,varargin)
%CROPMOVIE allows the user to crop the movie
%
% Syntax:
% croppedMovieData = cropMovie(movieData,cropRoi);
% croppedMovieData = cropMovie(movieData,cropRoi,newOutputDir,'additionalFiles',additionalFiles)
% 
% Description:
% 
% This function crops the channels of the input movie as well as any
% specified additional file according to the selected roi. The cropped
% channels are saved in a new directory and a new movie data is generated
% using the same settings as the original one losing any performed process
%
% Input: 
%
%   movieData - The MovieData object to be cropped
% 
%   cropROI - a four-element position vector [xmin ymin width height] that
%             specifies the size and position of the cropping rectangle.
% 
%   outputDirectory - (optional) A string containing the path where the
%   cropped movie should be saved. If not input, a dialog box will ask the
%   user to select a directory
%   
%   additionalFiles - (optional,param/value) a cell array of paths
%   corresponding to additional images to be cropped using the same region.

% Sebastien Besson, Sep 2011


% Input check
ip= inputParser;
ip.addRequired('movieData',@(x) isa(x,'MovieData'));
ip.addRequired('cropROI',@(x) isvector(x) && numel(x)==4);
ip.addOptional('outputDirectory','',@ischar)
ip.addParamValue('additionalFiles',{},@(x) iscell);
ip.parse(movieData,cropROI,varargin{:});
additionalFiles=ip.Results.additionalFiles;
if isempty(ip.Results.outputDirectory)
    outputDirectory = uigetdir(pwd,'Select output directory for cropped Movie');
    if isequal(outputDirectory,0), return; end
end

% Create log message
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing');
    logMsg = @(chan) ['Please wait, cropping images for channel ' num2str(chan)];
    timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
    tic;
end

% Read channel information (number of frames, channel names)
nFrames = movieData.nFrames_;
imDirs = movieData.getChannelPaths();
imageFileNames = movieData.getImageFileNames();
nChan = numel(movieData.channels_);
nTot = nChan*nFrames;
inImage = @(chan,frame) [imDirs{chan} filesep imageFileNames{chan}{frame}];

% Create new channel directory names for image writing
[~,chanNames]=cellfun(@fileparts,imDirs,'UniformOutput',false);
newImDirs = cellfun(@(x) [outputDirectory filesep x],chanNames,...
    'UniformOutput',false);
outImage = @(chan,frame) [newImDirs{chan} filesep imageFileNames{chan}{frame}];

% Read public access channel properties
m=?Channel;
channelFieldsAccess=cellfun(@(x) x.SetAccess,m.Properties,'Unif',false);
channelPublicFields= cellfun(@(x) strcmpi(x,'public'),channelFieldsAccess);

%Copy channel images
channels(nChan)=Channel();
for i = 1:nChan
    disp('Cropping channel:')
    disp(newImDirs{i});
    disp('Results will be saved under:')
    disp(newImDirs{i});
    mkClrDir(newImDirs{i});
    
    % Create channel object and copy public properties
    channels(i)=Channel(newImDirs{i});
    s= struct(movieData.channels_(i));
    fields=fieldnames(s);    
    set(channels(i),rmfield(s,fields(~channelPublicFields)));
    
    for j= 1:nFrames
        % Read original image, crop it and save it
        imwrite(imcrop(imread(inImage(i,j)),cropROI), outImage(i,j));
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(i) timeMsg(tj*nTot/nj-tj)]));
        end
    end
end

% Crop and write additional files at the base of the ouput directory
if ~isempty(additionalFiles)
    if ishandle(wtBar),
        waitbar(1,wtBar,'Please wait, cropping additional files...');
    end
    for i = 1:numel(additionalFiles)
        [~,fileName,fileExt]=fileparts(additionalFiles{i});

        % Read original image, crop it and save it
        imwrite(imcrop(imread(additionalFiles{i}),cropROI),...
            [outputDirectory filesep fileName fileExt]);
    end
end

if ishandle(wtBar),
    waitbar(1,wtBar,'Please wait, creating movieData object...');
end

% Read public access & unchanged movie properties
m=?MovieData;
movieFieldsAccess=cellfun(@(x) x.SetAccess,m.Properties,'Unif',false);
moviePublicFields= cellfun(@(x) strcmpi(x,'public'),movieFieldsAccess);
changedFields = {'outputDirectory_','movieDataPath_','movieDataFileName_'};
movieChangedFields= cellfun(@(x) any(strcmpi(x.Name,changedFields)),m.Properties);

% Create movieData object and copy public properties
croppedMovieData=MovieData(channels,outputDirectory,'movieDataPath_',outputDirectory,...
    'movieDataFileName_','movieData.mat');
s= struct(movieData);
fields=fieldnames(s);
set(croppedMovieData,rmfield(s,fields(~moviePublicFields | movieChangedFields)));

% Perform sanityCheck
croppedMovieData.sanityCheck

if ishandle(wtBar), close(wtBar); end
