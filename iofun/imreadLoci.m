function [image,metadata] = imreadLoci(fileName,options)
%IMREADLOCI reads images using the bio-formats java library
%
%SYNOPSIS: [image,metadata] = imreadLoci(fileName,options)
%
% INPUT fileName: Full filename including path of the image file you try to
%           open. If empty, imreadLoci will open a dialog
%       options : tbd
%
% OUTPUT image : image array
%        metadata : metadata - only .imageSize implemented
%
% REMARKS:  This function is not quite finished yet.
%           You need loci_tools.jar on your java path. 
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-Sep-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

image = [];
metadata = [];
        
%=========================
% TEST INPUT
%=========================

% check for fileName and load if necessary
if nargin == 0 || isempty(fileName)
    [fileName,pathName] = uigetfile({'*.*','All Files'},'Please select an image file');
    if fileName == 0
        % user aborted. Exit gracefully
        
        disp('No file selected')
        return
    else
        fileName = fullfile(pathName,fileName);
    end
else
    % check if file exists
    if ~exist(fileName,'file')
        error('file %s not found',fileName);
    end
end

% check for options
if nargin < 2 || isempty(options)
    options = [];
end

%============================



%============================
% LOAD IMAGE FILE
%============================

% find image name 
[imagePath, imageName] = fileparts(fileName);

% create reader object
reader = loci.formats.ImageReader;
% reader = loci.formats.ChannelFiller();
% reader = loci.formats.ChannelSeparator(reader);
% reader = loci.formats.FileStitcher(reader);

% supply the image
%progressText(0,sprintf('reading properties for %s',imageName))
reader.setId(fileName);

% read image size
imSize = zeros(1,5);
imSize(1) = reader.getSizeX;
imSize(2) = reader.getSizeY;
imSize(3) = reader.getSizeZ;
imSize(4) = reader.getSizeC;
imSize(5) = reader.getSizeT;

% -- insert size check here

% preassign arrays
image = zeros(imSize);
imagePlane = zeros(imSize(1:2));

progressText(1);

% maybe tiff-series are loaded this way?
% numSeries = r.getSeriesCount();
% for s = 1:numSeries
%     fprintf('Reading series #%d', s);
%     r.setSeries(s - 1);

% read images and store
progressText(0,sprintf('loading %s',imageName));

numImages = reader.getImageCount;
for i = 1:numImages
    img = reader.openImage(i - 1);
    % convert Java BufferedImage to MATLAB image
    imagePlane(:) = img.getData.getPixels(0, 0, imSize(1), imSize(2), []);

    % get current z/c/t
    zct = reader.getZCTCoords(i - 1) + 1;

    % store
    image(:,:,zct(1),zct(2),zct(3)) = imagePlane;
    
    progressText(i/numImages);

end

% read metadata
metadata.imSize = imSize;
