function MD = bfImport(dataPath,varargin)
% BFIMPORT imports movie files into MovieData objects using Bioformats 
%
% MD = bfimport(dataPath)
%
% Load proprietary files using the Bioformats library. Read the metadata
% that is associated with the movie and the channels and set them into the
% created movie objects. Optionally images can be extracted and saved as
% individual TIFF files.
%
% Input:
% 
%   dataPath - A string containing the full path to the movie file.
%
%   extractImages - Optional. If true, individual images will be extracted
%   and saved as TIF images.
%
% Output:
%
%   MD - A MovieData object

% Sebastien Besson, Dec 2011

status = bfCheckJavaPath();
assert(status, 'Bioformats library missing');

% Input check
ip=inputParser;
ip.addRequired('dataPath',@ischar);
ip.addOptional('extractImages',false,@islogical);
ip.addParamValue('outputDirectory',[],@ischar);
ip.parse(dataPath,varargin{:});
extractImages = ip.Results.extractImages;

assert(exist(dataPath,'file')==2,'File does not exist'); % Check path

try
    % Retrieve movie reader and metadata
    r=bfGetReader(dataPath);
    r.setSeries(0);
catch bfException
    ME = MException('lccb:import:error','Import error');
    ME = ME.addCause(bfException);
    throw(ME);
end

iSeries =0;
movieArgs = getMovieMetadata(r, iSeries);

% Read number of channels, frames and stacks
nFrames =  r.getMetadataStore().getPixelsSizeT(iSeries).getValue;
nChan =  r.getMetadataStore().getPixelsSizeC(iSeries).getValue;
nZ =  r.getMetadataStore().getPixelsSizeZ(iSeries).getValue;

% Set output directory (based on image extraction flag)
[mainPath,movieName]=fileparts(dataPath);
if ~isempty(ip.Results.outputDirectory)
    outputDir = ip.Results.outputDirectory;
else
    outputDir = fullfile(mainPath, movieName);
end
movieFileName=[movieName '.mat'];

% Create movie channels
channelPath=cell(1,nChan);
movieChannels(1,nChan)=Channel();
for iChan = 1:nChan

    channelArgs = getChannelMetadata(r, iSeries, iChan-1);
    % Read channelName
    chanName=r.getMetadataStore().getChannelName(iSeries, iChan-1);
    if isempty(chanName), 
        chanName = ['Channel_' num2str(iChan)]; 
    else
        chanName = char(chanName.toString); 
    end
    
    % Create new channel
    if extractImages
        channelPath{iChan} = fullfile(outputDir, chanName);
    else
        channelPath{iChan} = dataPath;
    end
    movieChannels(iChan) = Channel(channelPath{iChan}, channelArgs{:});
end

% Create movie object
MD=MovieData(movieChannels, outputDir, movieArgs{:});
MD.setPath(outputDir);
MD.setFilename(movieFileName);


if extractImages    
    % Create anonymous functions for reading files
    tString=@(t)num2str(t, ['%0' num2str(floor(log10(nFrames))+1) '.f']);
    zString=@(z)num2str(z, ['%0' num2str(floor(log10(nZ))+1) '.f']);
    imageName = @(c,t,z) [movieName '_w' num2str(movieChannels(c).emissionWavelength_) ...
        '_z' zString(z),'_t' tString(t),'.tif'];

    % Clean channel directories and save images as TIF files
    for iChan = 1:nChan, mkClrDir(channelPath{iChan}); end
    for iPlane = 1:r.getImageCount()
        index = r.getZCTCoords(iPlane - 1);
        imwrite(bfGetPlane(r, iPlane),[channelPath{index(2) + 1} filesep ...
            imageName(index(2) + 1, index(3) + 1, index(1) + 1)],'tif');
    end
end

% Close reader and check movie sanity
r.close;
MD.sanityCheck;

function movieArgs = getMovieMetadata(r, iSeries)

% Create movie metadata cell array using read metadata
movieArgs={};

pixelSizeX = r.getMetadataStore().getPixelsPhysicalSizeX(iSeries);
% Pixel size might be automatically set to 1.0 by @#$% Metamorph
hasValidPixelSize = ~isempty(pixelSizeX) && pixelSizeX.getValue ~= 1;
if hasValidPixelSize
    % Convert from microns to nm and check x and y values are equal
    pixelSizeX= pixelSizeX.getValue*10^3;
    pixelSizeY= r.getMetadataStore().getPixelsPhysicalSizeY(iSeries).getValue*10^3;
    assert(isequal(pixelSizeX,pixelSizeY),'Pixel size different in x and y');
    movieArgs=horzcat(movieArgs,'pixelSize_',pixelSizeX);
end

% Camera bit depth
camBitdepth = r.getBitsPerPixel();
hasValidCamBitDepth = ~isempty(camBitdepth) && mod(camBitdepth, 2) == 0;
if hasValidCamBitDepth
    movieArgs=horzcat(movieArgs,'camBitdepth_',camBitdepth);
end

% Time interval
timeInterval = r.getMetadataStore().getPixelsTimeIncrement(iSeries);
if ~isempty(timeInterval)
    movieArgs=horzcat(movieArgs,'timeInterval_',double(timeInterval));
end

% Lens numerical aperture
try % Use a try-catch statement because property is not always defined
    lensNA=r.getMetadataStore().getObjectiveLensNA(0,0);
    if ~isempty(lensNA)
        movieArgs=horzcat(movieArgs,'numAperture_',double(lensNA));
    elseif ~isempty(r.getMetadataStore().getObjectiveID(0,0))
        % Hard-coded for deltavision files. Try to get the objective id and
        % read the objective na from a lookup table
        tokens=regexp(char(r.getMetadataStore().getObjectiveID(0,0).toString),...
            '^Objective\:= (\d+)$','once','tokens');
        if ~isempty(tokens)
            [na,mag]=getLensProperties(str2double(tokens),{'na','magn'});
            movieArgs=horzcat(movieArgs,'numAperture_',na,'magnification_',mag);
        end
    end
end

function channelArgs = getChannelMetadata(r, iSeries, iChan)

channelArgs={};

% Read excitation wavelength
exwlgth=r.getMetadataStore().getChannelExcitationWavelength(iSeries, iChan);
if ~isempty(exwlgth)
    channelArgs=horzcat(channelArgs, 'excitationWavelength_', exwlgth.getValue);
end

% Fill emission wavelength
emwlgth=r.getMetadataStore().getChannelEmissionWavelength(iSeries, iChan);
if isempty(emwlgth)
    try
        emwlgth= r.getMetadataStore().getChannelLightSourceSettingsWavelength(iSeries, iChan);
    end
end
if ~isempty(emwlgth)
    channelArgs = horzcat(channelArgs, 'emissionWavelength_', emwlgth.getValue);
end