function MD = omeroCacheImport(session,imageID,varargin)
% OMEROCACHEIMPORT caches images from an OMERO server into MovieData objects
%
% movie = omeroCacheImport(session, imageID)
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
ip.addRequired('session',@MovieObject.isOmeroSession);
ip.addRequired('imageID',@isscalar);
ip.addOptional('outputDirectory',[],@ischar);
ip.parse(session,imageID,varargin{:});

if isempty(ip.Results.outputDirectory)
    [~, outputDir] = uiputfile('*.mat','Find a place to save your analysis',...
        'movieData.mat');
    if isequal(outputDir,0), return; end
else
    outputDir=ip.Results.outputDirectory;
    if ~isdir(outputDir), mkdir(outputDir); end
end

rawDataFile = fullfile(outputDir, [num2str(imageID) '.ome.tiff']);

exportImageAsOMETIFF(session, imageID, rawDataFile);

MD = MovieData.load(rawDataFile);
MD.setOmeroId(imageID);
MD.setOmeroSave(false);
MD.save;