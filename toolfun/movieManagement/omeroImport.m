function movie = omeroImport(image,session,varargin)
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
ip.addRequired('session',@(x) isa(x,'omero.api.ServiceFactoryPrxHelper'));
ip.addOptional('outputDirectory',[],@ischar);
ip.parse(image,session,varargin{:});

movieArgs={}; % Create properties cell array based on existing metadata

% Retrieve pixels
svc=session.getPixelsService();
pixels=svc.retrievePixDescription(image.getPixels(0).getId.getValue);

% Get pixel size
pixelSizeX = pixels.getPhysicalSizeX;
if ~isempty(pixelSizeX)
    pixelSizeX = pixelSizeX.getValue*10^3; % Convert to nm
    pixelSizeY = pixels.getPhysicalSizeY.getValue*10^3;
    assert(isequal(pixelSizeX,pixelSizeY),'Pixel size different in x and y');
    movieArgs=horzcat(movieArgs,'pixelSize_',pixelSizeX);
end

% Get camera bit depth
camBitdepth = pixels.getPixelsType.getBitSize.getValue/2;
if ~isempty(camBitdepth)
    movieArgs=horzcat(movieArgs,'camBitdepth_',camBitdepth);
end

% Get time interval
timeInterval = pixels.getTimeIncrement();
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
nFrames =  pixels.getSizeT.getValue;
nChan =  pixels.getSizeC.getValue;
nZ =  pixels.getSizeZ.getValue;

% Set output directory (based on image extraction flag)
movie=MovieData;

if isempty(ip.Results.outputDirectory)
    [movieFileName,outputDir] = uiputfile('*.mat','Find a place to save your analysis',...
        'movieData.mat');
    if isequal(outputDir,0), return; end
else
    outputDir=ip.Results.outputDirectory;
    if ~isdir(outputDir), mkdir(outputDir); end
    movieFileName='movie.mat';
end

% Create movie channels
movieChannels(1,nChan)=Channel();
channelArgs=cell(1,nChan);
channels = pixels.copyChannels;
for i=1:nChan
    channelArgs{i}={};

%     Read excitation wavelength
    exwlgth=channels.get(i-1).getLogicalChannel().getEmissionWave;
%     exwlgth=channels.get(i-1).getEmissionWave;

    if ~isempty(exwlgth) && exwlgth.getValue ~=1
        channelArgs{i}=horzcat(channelArgs{i},'excitationWavelength_',exwlgth.getValue);
    end
    
    % Fill emission wavelength
    emwlgth=channels.get(i-1).getLogicalChannel().getExcitationWave();

    if ~isempty(emwlgth) && emwlgth.getValue ~=1
        channelArgs{i}=horzcat(channelArgs{i},'emissionWavelength_',emwlgth.getValue);
    end
    
    % Read channelName
    chanName=channels.get(i-1).getLogicalChannel.getName;
    if isempty(chanName), 
        chanName = ['Channel_' num2str(i)]; 
    else
        chanName = char(chanName.toString); 
    end
    
    % Create new channel
%     channelPath{i}=dataPath;
%     end
    movieChannels(i)=Channel(pixels,channelArgs{i}{:});
end

% Create movie object
movie=MovieData(movieChannels,outputDir,movieArgs{:});
movie.setPath(outputDir);
movie.setFilename(movieFileName);
movie.setSession(session);

movie.sanityCheck;