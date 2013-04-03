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
ip.addRequired('session',@(x) isa(x,'omero.api.ServiceFactoryPrxHelper'));
ip.addRequired('imageID',@isscalar);
ip.addOptional('outputDirectory',[],@ischar);
ip.parse(session,imageID,varargin{:});

% Get images
proxy = session.getContainerService();
id = toJavaList(imageID, 'java.lang.Long');
imageList = proxy.getImages('omero.model.Image', id, omero.sys.ParametersI());
images = toMatlabList(imageList);
assert(~isempty(images) && ~isempty(images(1)), 'No image found');
image=images(1);

% Retrieve pixels
svc=session.getPixelsService();
pixels=svc.retrievePixDescription(image.getPixels(0).getId.getValue);


movieArgs={}; % Create properties cell array based on existing metadata



pixelsA=pojos.PixelsData(pixels);

% Get pixel size
pixelSize = pixelsA.getPixelSizeX;
if ~isempty(pixelSize)
    assert(isequal(pixelSize,pixelsA.getPixelSizeY),'Pixel size different in x and y');
%     movieArgs=horzcat(movieArgs,'pixelSize_',pixelSize*1e3);
end

% Get camera bit depth
camBitdepth = pixels.getPixelsType.getBitSize.getValue/2;
if ~isempty(camBitdepth)
    movieArgs=horzcat(movieArgs,'camBitdepth_',camBitdepth);
end

% Get time interval
timeInterval = pixels.getTimeIncrement();
if ~isempty(timeInterval)
    movieArgs=horzcat(movieArgs,'timeInterval_', timeInterval.getValue());
end

% Get the lens numerical aperture
imageA = pojos.ImageAcquisitionData(image);
objective = imageA.getObjective;
if ~isempty(objective)
    lensNA = imageA.getObjective.getLensNA;
    movieArgs=horzcat(movieArgs,'numAperture_',double(lensNA));
end


% Read number of channels, frames and stacks
nChan =  pixelsA.getSizeC;

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

    channelD = pojos.ChannelData(i-1,channels.get(i-1));
    % Read excitation wavelength
    exwlgth=channelD.getExcitationWavelength();
    if ~isempty(exwlgth) && exwlgth ~= -1
        channelArgs{i}=horzcat(channelArgs{i},'excitationWavelength_',exwlgth);
    end
    
    % Read emission wavelength
    emwlgth=channelD.getEmissionWavelength();
    if ~isempty(emwlgth) && emwlgth ~= -1 && emwlgth ~= 1 % Bug
        channelArgs{i}=horzcat(channelArgs{i},'emissionWavelength_',emwlgth);
    end
    
    % Read channelName
    chanName=channelD.getName;
    if isempty(chanName), chanName = ['Channel_' num2str(i)]; end
    
    % Create new channel
%     channelPath{i}=dataPath;
%     end
    movieChannels(i)=Channel('',channelArgs{i}{:});
end

% Create movie object
movie=MovieData(movieChannels,outputDir,movieArgs{:});
movie.setPath(outputDir);
movie.setFilename(movieFileName);
movie.setOmeroId(imageID);
movie.setSession(session);
movie.setOmeroSave(true);

movie.sanityCheck;