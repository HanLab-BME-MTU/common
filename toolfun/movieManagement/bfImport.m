function movie = bfImport(dataPath,varargin)
% Import properiary movie files and convert them into MovieData objects
%
% movie = bfimport(dataPath)
%
% Load proprietary files using the bioformats library. Read the metadata
% that is associated with the movie and the channels and set them into the
% created movie objects. Optionally extract images.
%
% Input:
% 
%   dataPath - A string containing the full path to the movie proprietary file.
%
%   extractImages - Optional. A boolean to extract movie images and set
%   them into 
%
% Output:
%
%   movie - A MovieData object

% Sebastien Besson, Dec 2011

assert(exist(which('bfopen'),'file')==2,'Bioformats library missing');

% Input check
ip=inputParser;
ip.addRequired('dataPath',@ischar);
ip.addOptional('extractImages',true,@islogical)
ip.parse(dataPath,varargin{:});
extractImages = ip.Results.extractImages;

assert(exist(dataPath,'file')==2,'File does not exist');

try
    r=bfGetReader(dataPath);
    metadata=r.getMetadataStore();
catch ME
    ME2 = MException('lccb:import:error','Import error');
    ME2.addCause(ME);
    throw(ME2);
end



movieArgs={};

% Read pixel size
pixelSize = metadata.getPixelsPhysicalSizeX(0);
if ~isempty(pixelSize)
    assert(isequal(pixelSize.getValue,metadata.getPixelsPhysicalSizeY(0).getValue));
    movieArgs=horzcat(movieArgs,'pixelSize_',pixelSize.getValue*10^3);
end

% Read camera bit depath size
camBitdepth = r.getBitsPerPixel;
if ~isempty(camBitdepth)
    movieArgs=horzcat(movieArgs,'camBitdepth_',camBitdepth);
end


% Read time interval
timeInterval = metadata.getPixelsTimeIncrement(0);
if ~isempty(timeInterval)
    movieArgs=horzcat(movieArgs,'timeInterval_',double(timeInterval));
end

% Read the lens numerical aperture
try % Ue a tr-catch statement because property is not always defined
    lensNA=metadata.getObjectiveLensNA(0,0);
    if ~isempty(lensNA)
        movieArgs=horzcat(movieArgs,'numAperture_',double(lensNA));
    elseif ~isempty(metadata.getObjectiveID(0,0))
        % Hard-coded for deltavision files. Try to get the objective Id and
        % read the na from a lookup table
        tokens=regexp(char(metadata.getObjectiveID(0,0).toString),...
            '^Objective\:= (\d+)$','once','tokens');
        if ~isempty(tokens)
            movieArgs=horzcat(movieArgs,'numAperture_',naFromLensID(str2double(tokens)));
        end
    end
end

% Read number of channels, frames and stacks
[mainPath,movieName]=fileparts(dataPath);
if extractImages
    outputDir=[mainPath filesep movieName];
else
    outputDir=mainPath;
end
nFrames =  metadata.getPixelsSizeT(0).getValue;
nChan =  metadata.getPixelsSizeC(0).getValue;

% Create channel objects
channelPath=cell(nChan,1);
movieChannels(nChan,1)=Channel();
channelArgs=cell(nChan,1);
for i=1:nChan
    channelArgs{i}={};

    % Read excitation wavelength
    exwlgth=metadata.getChannelExcitationWavelength(0,i-1);
    if ~isempty(exwlgth)
        channelArgs{i}=horzcat(channelArgs{i},'excitationWavelength_',exwlgth.getValue);
    end
    
    % Fill emission wavelength
    emwlgth=metadata.getChannelEmissionWavelength(0,i-1);
    if isempty(emwlgth)
        try
            emwlgth= metadata.getChannelLightSourceSettingsWavelength(0,0);
        end
    end
    if ~isempty(emwlgth)
        channelArgs{i}=horzcat(channelArgs{i},'emissionWavelength_',emwlgth.getValue);
    end
    
    % Read channelName
    chanName=metadata.getChannelName(0,i-1);
    if isempty(chanName), 
        chanName = ['Channel_' num2str(i)]; 
    else
        chanName = char(chanName.toString); 
    end
    
    % Create new channel
    if extractImages
        channelPath{i} = [outputDir filesep chanName];
    else
        channelPath{i}=dataPath;
    end
    movieChannels(i)=Channel(channelPath{i},channelArgs{i}{:});
end

% Create movie object
movie=MovieData(movieChannels,outputDir,movieArgs{:});
movie.setPath(outputDir);
movie.setFilename([movieName '.mat']);

% Save images
dimensionOrder =char(metadata.getPixelsDimensionOrder(0));
dimensions = arrayfun(@(x) metadata.(['getPixelsSize' x])(0).getValue,...
    dimensionOrder(3:end));

% Create anonymous functions for reading files
chanIndex = @(index) index(dimensionOrder(3:end)=='C');
zIndex = @(index) index(dimensionOrder(3:end)=='Z');
tIndex = @(index) index(dimensionOrder(3:end)=='T');
tString=@(t)num2str(t, ['%0' num2str(floor(log10(nFrames))+1) '.f']);
imageName = @(c,t) [movieName '-w' num2str(movieChannels(c).emissionWavelength_) ...
    '_t' tString(t),'.tif'];

% Save images
if extractImages
    for i=1:nChan, mkClrDir(channelPath{i}); end
    for iPlane = 1:r.getImageCount()
        [index(1),index(2),index(3)]=ind2sub(dimensions,iPlane);
        
        imwrite(bfGetPlane(r,iPlane),[channelPath{chanIndex(index)} filesep ...
            imageName(chanIndex(index),tIndex(index))],'tif');
    end
end

r.close;
movie.sanityCheck;
