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
userId = movieData.getSession().getAdminService().getEventContext().userId;
nsToInclude = java.util.ArrayList;
nsToInclude.add(namespace);
nsToExclude = java.util.ArrayList;

options = omero.sys.ParametersI;
options.exp(omero.rtypes.rlong(userId)); %load the annotation for a given user.
options.exp(omero.rtypes.rlong(userId)); %load the annotation for a given user.
metadataService = movieData.getSession().getMetadataService();
% retrieve the annotations linked to images, for datasets use: 'omero.model.Dataset'
annotations = metadataService.loadSpecifiedAnnotations('omero.model.FileAnnotation', nsToInclude, nsToExclude, options);

% imageIds = java.util.ArrayList;
% imageIds.add(java.lang.Long(MD.getImage().getId().getValue))
% annotations = metadataService.loadAnnotation('omero.model.Image', imageIds,...
%     java.util.ArrayList,java.util.ArrayList,omero.sys.ParametersI());

if annotations.size() >0
    originalFile = annotations.get(0).getFile();
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

if annotations.size()==0
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