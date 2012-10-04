function omeroExportDetection(movieData,movieInfo)
% omeroExportDetection export the output of a detection to the OMERO server
%
% omeroExportDetection export the output of any detection process
% (generating a movieInfo like structure into OMERO. After creates a single
% ROI, the function creates as many Point shapes as the number of detected
% objects and sets the timepoints and x,y coordinates of this point shape
% using the detection information. Finally, the roi is attached to the
% image on the server.
%
% omeroExportDetection(movieData,movieInfo)
%
% Input:
% 
%   movieData - A MovieData object
%
%   movieInfo - The output of a detection process
%
% Output:
%

% Sebastien Besson, Jun 2012 (last modified Oct 2012)

ns = 'hms-detection';

% Input check
ip=inputParser;
ip.addRequired('movieData',@(x) isa(x,'MovieData') && x.isOmero());
ip.addRequired('movieInfo',@isstruct);
ip.parse(movieData,movieInfo);

% Retrieve image and update service
image = movieData.getImage();
roiService = movieData.getSession().getRoiService();
updateService = movieData.getSession().getUpdateService();

% Get previously saved ROIs with same namespace
roiOptions = omero.api.RoiOptions();
roiOptions.namespace = omero.rtypes.rstring(ns);
rois = roiService.findByImage(movieData.omeroId_, roiOptions).rois();

if rois.size()> 0
    % List rois to remove
    fprintf('Deleting %g existing rois with namespace %s\n', rois.size(), ns);
    list = javaArray('omero.api.delete.DeleteCommand', 1);
    for i = 1:rois().size()
        roiId = rois.get(i-1).getId().getValue;
        list(i) = omero.api.delete.DeleteCommand('/Roi', roiId, []);
    end
    
    %Delete the rois
    movieData.getSession().getDeleteService().queueDelete(list);
end

% Create ROI to attach to the image
progressText(0, 'Exporting ROIs')
roi = omero.model.RoiI();
roi.setImage(image);
roi.setNamespaces(ns);

for t=1:size(movieInfo,1)
    
    for i = 1:size(movieInfo(t).xCoord,1)        
        % Create point shape
        point = omero.model.PointI;
        point.setCx(omero.rtypes.rdouble(movieInfo(t).xCoord(i,1)));
        point.setCy(omero.rtypes.rdouble(movieInfo(t).yCoord(i,1)));
        point.setTheT(omero.rtypes.rint(t-1));
        point.setTheZ(omero.rtypes.rint(0));
        roi.addShape(point);
    end

    progressText(t/size(movieInfo,1))
end
updateService.saveAndReturnObject(roi);