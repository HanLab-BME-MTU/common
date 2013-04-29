function omeroSave(movieObject)
% OMEROSAVE uploads the output directory into OMERO as a file annotation
%
% omeroSave first create a zipped archive of all the content of the movie
% output directory. It then looks for a file annotation with the correct
% namespace attached to the Image. If existing, it uses the corresponding
% file else it create a new OriginalFile, saves the content of the zip file
% into this OriginalFile, uploads it to the server. Finally if a new file
% has been created, a new file annotation linking it to the image is
% created and uploaded to the server.
%
% omeroSave(movieObject)
%
% Input:
% 
%   movieObject - A MovieData object
%

% Sebastien Besson, Jun 2012 (last modified Oct 2012)

% To be replaced by omero.constants....
namespace = 'hms-tracking';
zipName = 'HMS-tracking.zip';

% Input check
ip=inputParser;
ip.addRequired('movieObject',@(x) isa(x,'MovieData') && x.isOmero() && x.canUpload());
ip.parse(movieObject);

% Zip output directory for attachment
zipPath = fileparts(movieObject.outputDirectory_);
zipFullPath = fullfile(zipPath,zipName);
zip(zipFullPath, movieObject.outputDirectory_)

% Create java io File
file = java.io.File(zipFullPath);
name = file.getName();
absolutePath = file.getAbsolutePath();
path = absolutePath.substring(0, absolutePath.length()-name.length());

% Load existing file annotations
fas = getImageFileAnnotations(movieObject.getSession(), movieObject.omeroId_,...
    'include', namespace);

if ~isempty(fas)
    % Read file of first found file annotation
    originalFile = fas(1).getFile;
else
    originalFile = omero.model.OriginalFileI;
end
originalFile.setName(rstring(name));
originalFile.setPath(rstring(path));
originalFile.setSize(rlong(file.length()));
originalFile.setSha1(rstring(''));
originalFile.setMimetype(rstring('application/zip'));

% now we save the originalFile object
iUpdate = movieObject.getSession().getUpdateService();
originalFile = iUpdate.saveAndReturnObject(originalFile);

% Initialize the service to load the raw data
rawFileStore = movieObject.getSession().createRawFileStore();
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

if isempty(fas)
    % now we have an original File in DB and raw data uploaded.
    % We now need to link the Original file to the image using the File annotation object. That's the way to do it.
    fa = omero.model.FileAnnotationI;
    fa.setFile(originalFile);
    fa.setDescription(rstring('HMS tracking')); % The description set above e.g. PointsModel
    fa.setNs(rstring(namespace)) % The name space you have set to identify the file annotation.
    
    % save the file annotation.
    fa = iUpdate.saveAndReturnObject(fa);
    
    % now link the image and the annotation
    link = omero.model.ImageAnnotationLinkI;
    link.setChild(fa);
    images = getImages(movieObject.getSession(), movieObject.omeroId_);
    if ~isempty(images)
        link.setParent(images(1));
        % save the link back to the server.
        iUpdate.saveAndReturnObject(link);
    end
end

% Delete zip file
delete(zipFullPath);