function movie = omeroImport(session,imageID,varargin)
% OMEROIMPORT imports images from an OMERO server into MovieData objects
%
% movie = omeroImport(session)
%
% Load proprietary files using the Bioformats library. Read the metadata
% that is associated with the movie and the channels and set them into the
% created movie objects. Optionally images can be extracted and saved as
% individual TIFF files.
%
% Input:
% 
%   session - an omero session
%
%   imageID - A string containing the full path to the movie file.
%
%   extractImages - Optional. If true, individual images will be extracted
%   and saved as TIF images.
%
% Output:
%
%   movie - A MovieData object

% Sebastien Besson, Dec 2011 (last modified Nov 2012)

% Input check
ip=inputParser;
ip.addRequired('session', @(x) isa(x,'omero.api.ServiceFactoryPrxHelper'));
ip.addRequired('imageID', @isscalar);
ip.addOptional('outputDirectory',[],@ischar);
ip.parse(session,imageID,varargin{:});

% Retrieve image and pixels
image = getImages(session, imageID);
assert(~isempty(image), 'No image of id %g found', imageID);
pixels = image.getPrimaryPixels();

% Create properties cell array based on existing metadata
metadataService = session.getMetadataService;
movieArgs={}; 

% Read pixel size
pixelSize = pixels.getPhysicalSizeX;
if ~isempty(pixelSize)
    assert(isequal(pixelSize, pixels.getPhysicalSizeY),...
        'Pixel size different in x and y');
    movieArgs = horzcat(movieArgs, 'pixelSize_', pixelSize.getValue * 1e3);
end

% Read camera bit depth
camBitdepth = pixels.getPixelsType.getBitSize.getValue/2;
if ~isempty(camBitdepth)
    movieArgs=horzcat(movieArgs,'camBitdepth_',camBitdepth);
end

% Read time interval
timeInterval = pixels.getTimeIncrement();
if ~isempty(timeInterval)
    movieArgs=horzcat(movieArgs,'timeInterval_', timeInterval.getValue());
end

% Read the lens numerical aperture
objectiveSettings = image.getObjectiveSettings;
if ~isempty(objectiveSettings)
    instrument = metadataService.loadInstrument(objectiveSettings.getId.getValue);
    objective = instrument.copyObjective().get(0);
    if ~isempty(objective)
        lensNA = objective.getLensNA().getValue();
        movieArgs=horzcat(movieArgs,'numAperture_',double(lensNA));
    end
end


% Read number of channels, frames and stacks
nChan =  pixels.getSizeC().getValue();

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
pixelsId = toJavaList(pixels.getId.getValue,'java.lang.Long');
omeroChannels = metadataService.loadChannelAcquisitionData(pixelsId);
for i=1:nChan
    channelArgs{i}={};
    omeroChannel = omeroChannels.get(i-1);
    
    % Read excitation wavelength
    exwlgth = omeroChannel.getExcitationWave().getValue();
    if ~isempty(exwlgth) && exwlgth ~= -1
        channelArgs{i}=horzcat(channelArgs{i},'excitationWavelength_',exwlgth);
    end
    
    % Read emission wavelength
    emwlgth=omeroChannel.getEmissionWave().getValue();
    if ~isempty(emwlgth) && emwlgth ~= -1 && emwlgth ~= 1 % Bug
        channelArgs{i}=horzcat(channelArgs{i},'emissionWavelength_',emwlgth);
    end
    
    % Read channel xame
    chanName = omeroChannels.get(0).getName.getValue;
    if isempty(chanName), chanName = ['Channel_' num2str(i)]; end
    
    movieChannels(i)=Channel('',channelArgs{i}{:});
end

% Create movie object
movie=MovieData(movieChannels,outputDir,movieArgs{:});
movie.setPath(outputDir);
movie.setFilename(movieFileName);
movie.setOmeroId(imageID);
movie.setOmeroSession(session);
movie.setOmeroSave(true);

movie.sanityCheck;