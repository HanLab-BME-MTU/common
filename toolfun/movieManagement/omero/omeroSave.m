function omeroSave(movieData)
% omeroSave uploads the outputDirectory into OMERO as a file annotation
%
% omeroSave first create a zipped archive of all the content of the movie
% output directory. It then looks for a file annotation with the correct
% namespace attached to the Image. If existing, it uses the corresponding
% file else it create a new OriginalFile, saves the content of the zip file
% into this OriginalFile, uploads it to the server. Finally if a new file
% has been created, a new file annotation linking it to the image is
% created and uploaded to the server.
%
% omeroSave(movieData)
%
% Input:
% 
%   movieData - A MovieData object
%

% Sebastien Besson, Jun 2012 (last modified Oct 2012)

% To be replaced by omero.constants....
namespace = 'hms-tracking';
zipName = 'HMS-tracking.zip';

% Input check
ip=inputParser;
ip.addRequired('movieData',@(x) isa(x,'MovieData') && x.isOmero());
ip.parse(movieData);

% Zip output directory for attachment
zipPath = fileparts(movieData.outputDirectory_);
zipFullPath = fullfile(zipPath,zipName);
zip(zipFullPath, movieData.outputDirectory_)

% Create java io File
file = java.io.File(zipFullPath);
name = file.getName();
absolutePath = file.getAbsolutePath();
path = absolutePath.substring(0, absolutePath.length()-name.length());

% Load existing file annotations
files = getOmeroFileAnnotations(movieData.getSession(), movieData.omeroId_);

if ~isempty(files{1})
    originalFile = files{1};
else
    originalFile = omero.model.OriginalFileI;
end
originalFile.setName(omero.rtypes.rstring(name));
originalFile.setPath(omero.rtypes.rstring(path));
originalFile.setSize(omero.rtypes.rlong(file.length()));
originalFile.setSha1(omero.rtypes.rstring('a'));
originalFile.setMimetype(omero.rtypes.rstring('application/zip'));

% now we save the originalFile object
iUpdate = movieData.getSession().getUpdateService();
originalFile = iUpdate.saveAndReturnObject(originalFile);

% Initialize the service to load the raw data
rawFileStore = movieData.getSession().createRawFileStore();
rawFileStore.setFileId(originalFile.getId().getValue());

%  open file and read it

%code for small file.
fid = fopen(zipFullPath);
byteArray = fread(fid,[1, file.length()], 'uint8');
rawFileStore.write(byteArray, 0, file.length());
fclose(fid);

originalFile = rawFileStore.save();
% Important to close the service
rawFileStore.close();

if isempty(files{1})
    % now we have an original File in DB and raw data uploaded.
    % We now need to link the Original file to the image using the File annotation object. That's the way to do it.
    fa = omero.model.FileAnnotationI;
    fa.setFile(originalFile);
    fa.setDescription(omero.rtypes.rstring('HMS tracking')); % The description set above e.g. PointsModel
    fa.setNs(omero.rtypes.rstring(namespace)) % The name space you have set to identify the file annotation.
    
    % save the file annotation.
    fa = iUpdate.saveAndReturnObject(fa);
    
    % now link the image and the annotation
    link = omero.model.ImageAnnotationLinkI;
    link.setChild(fa);
    link.setParent(movieData.getImage());
    % save the link back to the server.
    iUpdate.saveAndReturnObject(link);
end

% Delete zip file
delete(zipFullPath);