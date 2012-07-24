function omeroSave(movieData)
% omeroSave export the output of a detection into the OMERO model
%
% omeroSave(movieData,movieInfo)
%
% Input:
% 
%   movieData - A MovieObject object
%
% Output:
%

% Sebastien Besson, Jun 2012


% Input check
ip=inputParser;
ip.addRequired('movieData',@(x) isa(x,'MovieData'));
ip.parse(movieData);

image = movieData.getImage();
iUpdate = movieData.getSession().getUpdateService();

zipPath = movieData.outputDirectory_;
zipName = 'HMS-tracking.zip';
zipFullPath = fullfile(zipPath,zipName);
zip(fullfile(zipPath,zipName),movieData.outputDirectory_)

file = java.io.File(zipFullPath);
name = file.getName();
absolutePath = file.getAbsolutePath();
path = absolutePath.substring(0, absolutePath.length()-name.length());

originalFile = omero.model.OriginalFileI;
originalFile.setName(omero.rtypes.rstring(name));
originalFile.setPath(omero.rtypes.rstring(path));
originalFile.setSize(omero.rtypes.rlong(file.length()));
originalFile.setSha1(omero.rtypes.rstring('a'));
originalFile.setMimetype(omero.rtypes.rstring('application/zip'));

% now we save the originalFile object
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
% now we have an original File in DB and raw data uploaded.
% We now need to link the Original file to the image using the File annotation object. That's the way to do it.
fa = omero.model.FileAnnotationI;
fa.setFile(originalFile);
fa.setDescription(omero.rtypes.rstring('HMS tracking')); % The description set above e.g. PointsModel
fa.setNs(omero.rtypes.rstring('hms-tracking')) % The name space you have set to identify the file annotation.

% save the file annotation.
fa = iUpdate.saveAndReturnObject(fa);

% now link the image and the annotation
link = omero.model.ImageAnnotationLinkI;
link.setChild(fa);
link.setParent(image);
% save the link back to the server.
iUpdate.saveAndReturnObject(link);

delete(zipFullPath);