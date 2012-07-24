function movie = omeroImport(image,varargin)
% BFIMPORT imports movie files into MovieData objects using Bioformats 
%
% movie = bfimport(image)
%
% Load proprietary files using the Bioformats library. Read the metadata
% that is associated with the movie and the channels and set them into the
% created movie objects. Optionally images can be extracted and saved as
% individual TIFF files.
%
% Input:
% 
%   image - A string containing the full path to the movie file.
%
%   extractImages - Optional. If true, individual images will be extracted
%   and saved as TIF images.
%
% Output:
%
%   movie - A MovieData object

% Sebastien Besson, Dec 2011

if ~exist('omero.client','class'), loadOmero; end

% Input check
ip=inputParser;
ip.addRequired('image',@(x) isa(x,'omero.model.ImageI'));
ip.addParamValue('outputDirectory',[],@ischar);
ip.parse(image,varargin{:});

movieArgs={}; % Create properties cell array based on existing metadata

% Get pixel size
pixelSizeX = image.getPixels(0).getPhysicalSizeX;
if ~isempty(pixelSizeX)
    pixelSizeX = pixelSizeX.getValue*10^3; % Convert to nm
    pixelSizeY = image.getPixels(0).getPhysicalSizeY.getValue*10^3;
    assert(isequal(pixelSizeX,pixelSizeY),'Pixel size different in x and y');
    movieArgs=horzcat(movieArgs,'pixelSize_',pixelSizeX);
end

% Get camera bit depth
camBitdepth = image.getPixels(0).getPixelsType.getBitSize.getValue/2;
if ~isempty(camBitdepth)
    movieArgs=horzcat(movieArgs,'camBitdepth_',camBitdepth);
end

% Get time interval
timeInterval = image.getPixels(0).getTimeIncrement();
if ~isempty(timeInterval)
    movieArgs=horzcat(movieArgs,'timeInterval_',double(timeInterval));
end

% Get the lens numerical aperture
try % Use a tyr-catch statement because property is not always defined
    lensNA=metadata.getObjectiveLensNA(0,0);
    if ~isempty(lensNA)
        movieArgs=horzcat(movieArgs,'numAperture_',double(lensNA));
    elseif ~isempty(metadata.getObjectiveID(0,0))
        % Hard-coded for deltavision files. Try to get the objective id and
        % read the objective na from a lookup table
        tokens=regexp(char(metadata.getObjectiveID(0,0).toString),...
            '^Objective\:= (\d+)$','once','tokens');
        if ~isempty(tokens)
            [na,mag]=getLensProperties(str2double(tokens),{'na','magn'});
            movieArgs=horzcat(movieArgs,'numAperture_',na,'magnification_',mag);
        end
    end
end

% Read number of channels, frames and stacks
nFrames =  image.getPixels(0).getSizeT.getValue;
nChan =  image.getPixels(0).getSizeC.getValue;
nZ =  image.getPixels(0).getSizeZ.getValue;

% Set output directory (based on image extraction flag)
movie=MovieData;

if isempty(ip.Results.outputDirectory)
    [movieFileName,outputDir] = uiputfile('*.mat','Find a place to save your analysis',...
        'movieData.mat');
    if isequal(outputDir,0), return; end
else
    outputDir=ip.Results.outputDirectory;
    if ~isdir(outputDir), mkdir(outputDir); end
    movieFileName=[movieName '.mat'];
end

% Create movie channels
channelPath=cell(1,nChan);
movieChannels(1,nChan)=Channel();
channelArgs=cell(1,nChan);
for i=1:nChan
    channelArgs{i}={};

    % Read excitation wavelength
%     exwlgth=metadata.getChannelExcitationWavelength(0,i-1);
%     if ~isempty(exwlgth)
%         channelArgs{i}=horzcat(channelArgs{i},'excitationWavelength_',exwlgth.getValue);
%     end
%     
%     % Fill emission wavelength
%     emwlgth=metadata.getChannelEmissionWavelength(0,i-1);
%     if isempty(emwlgth)
%         try
%             emwlgth= metadata.getChannelLightSourceSettingsWavelength(0,0);
%         end
%     end
%     if ~isempty(emwlgth)
%         channelArgs{i}=horzcat(channelArgs{i},'emissionWavelength_',emwlgth.getValue);
%     end
    
    % Read channelName
%     chanName=metadata.getChannelName(0,i-1);
%     if isempty(chanName), 
%         chanName = ['Channel_' num2str(i)]; 
%     else
%         chanName = char(chanName.toString); 
%     end
    
    % Create new channel
%     channelPath{i}=dataPath;
%     end
    movieChannels(i)=Channel(image.getPixels(0),channelArgs{i}{:});
end

% Create movie object
movie=MovieData(movieChannels,outputDir,movieArgs{:});
movie.setPath(outputDir);
movie.setFilename(movieFileName);

% 
% if extractImages    
%     % Get dimensions
%     dimensionOrder =char(metadata.getPixelsDimensionOrder(0));
%     dimensions = arrayfun(@(x) metadata.(['getPixelsSize' x])(0).getValue,...
%         dimensionOrder(3:end));
%     
%     % Create anonymous functions for reading files
%     chanIndex = @(index) index(dimensionOrder(3:end)=='C');
%     zIndex = @(index) index(dimensionOrder(3:end)=='Z');
%     tIndex = @(index) index(dimensionOrder(3:end)=='T');
%     tString=@(t)num2str(t, ['%0' num2str(floor(log10(nFrames))+1) '.f']);
%     zString=@(z)num2str(z, ['%0' num2str(floor(log10(nZ))+1) '.f']);
%     imageName = @(c,t,z) [movieName '_w' num2str(movieChannels(c).emissionWavelength_) ...
%         '_z' zString(z),'_t' tString(t),'.tif'];
% 
%     % Clean channel directories and save images as TIF files
%     for i=1:nChan, mkClrDir(channelPath{i}); end
%     for iPlane = 1:r.getImageCount()
%         [index(1),index(2),index(3)]=ind2sub(dimensions,iPlane);
%         
%         imwrite(bfGetPlane(r,iPlane),[channelPath{chanIndex(index)} filesep ...
%             imageName(chanIndex(index),tIndex(index),zIndex(index))],'tif');
%     end
% end

% Close reader and check movie sanity
% r.close;
movie.sanityCheck;