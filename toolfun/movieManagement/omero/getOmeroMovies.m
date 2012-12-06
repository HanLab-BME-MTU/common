function MD = getOmeroMovies(client, imageIDs, varargin)
% getOmeroMovies creates or loads MovieData object from OMERO images
%
% SYNOPSIS
%      
%
% INPUT
%    client -  an omero.client object with a created session
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
% Sebastien Besson, Nov 2012

% Input check
ip = inputParser;
ip.addRequired('imageIDs', @(x) isvector(x) || isa(x,'java.util.ArrayList'));
ip.addOptional('path', fullfile(getenv('HOME'), 'omero'), @ischar);
ip.parse(imageIDs, varargin{:});

% Create a java array list for the IDs
if ~isa(imageIDs, 'java.util.ArrayList')
    ids = java.util.ArrayList();
    for i = imageIDs(:)'
        ids.add(java.lang.Long(i)); %add the id of the image.
    end
    imageIDs = ids;
end

% Initialize movie array
nMovies = imageIDs.size;
MD(nMovies) = MovieData();

% Retrieve existing file annotations with the correct namespace
fileAnnotations = getOmeroFileAnnotations(client.getSession(), imageIDs);
hasFileAnnotation = ~cellfun(@isempty, fileAnnotations);

%% OMERO Images for which a MovieData has been created and uploaded
if any(hasFileAnnotation)
    % Set temporary file to extract file annotations
    zipPath = fullfile(ip.Results.path, 'tmp.zip');
    
    % List found image IDs and corresponding OriginalFile IDs
    foundIDs = find(hasFileAnnotation);
    fileIDs = cellfun(@(x) x.getId.getValue, fileAnnotations(hasFileAnnotation));
    
    for i = 1:numel(fileIDs)        
        iMovie = (foundIDs(i));
        
        % Download file annotation using omero.client.download
        localfile = java.io.File(zipPath);
        client.download(fileIDs(i), localfile);
        
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
            MD(iMovie) = MovieData.load(matFiles{j}, client.getSession, false);
        end
        
        if isempty(MD(iMovie)), error('No movie found'); end     
    end
end

%% OMERO Images where a new MovieData object needs to be created
if ~all(hasFileAnnotation)
    newIDs = find(~hasFileAnnotation);
    
    %% Retrieve a given plane.
    for i = newIDs
        path = fullfile(ip.Results.path, num2str(imageIDs.get(i-1)));        
        MD(i) = omeroImport(client.getSession(),imageIDs.get(i-1),path);
    end
end
