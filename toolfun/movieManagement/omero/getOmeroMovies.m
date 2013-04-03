function MD = getOmeroMovies(session, imageIDs, varargin)
% getOmeroMovies creates or loads MovieData object from OMERO images
%
% SYNOPSIS
%
%
% INPUT
%    session -  a session
%
%    imageIDs - an array of imageIDs. May be a Matlab array or a Java
%    ArrayList.
%
%    path - Optional. The default path where to extract/create the
%    MovieData objects for analysis
%
% OUTPUT
%    MD - an array of MovieData object corresponding to the images.
%
% Sebastien Besson, Nov 2012 (last modified Mar 2013)

% Input check
ip = inputParser;
ip.addRequired('imageIDs', @isvector);
ip.addOptional('path', fullfile(getenv('HOME'), 'omero'), @ischar);
ip.parse(imageIDs, varargin{:});

% Initialize movie array
nMovies = numel(imageIDs);
MD(nMovies) = MovieData();

% Set temporary file to extract file annotations
zipPath = fullfile(ip.Results.path, 'tmp.zip');

for i = 1 : nMovies
    newID = imageIDs(i);
    fas = getOmeroFileAnnotations(session, newID);
    
    if isempty(fas)
        path = fullfile(ip.Results.path, num2str(newID));
        MD(i) = omeroImport(session, newID, path);
    else
        getFileAnnotationContent(session, fas(1), zipPath);
        
        % Unzip and delete temporary fil
        zipFiles = unzip(zipPath, ip.Results.path);
        delete(zipPath);
        
        % List unzipped MAT files
        isMatFile = cellfun(@(x) strcmp(x(end-2:end),'mat'), zipFiles);
        matFiles = zipFiles(isMatFile);
        for j = 1: numel(matFiles)
            % Find MAT file containing MovieData object
            vars = whos('-file', matFiles{j});
            hasMovie = any(cellfun(@(x) strcmp(x, 'MovieData'),{vars.class}));
            if ~hasMovie, continue; end
            
            % Load MovieData object
            MD(i) = MovieData.load(matFiles{j}, session, false);
        end
    end
end

