function omeroExportDetection(movieData,movieInfo)
% omeroExportDetection export the output of a detection into the OMERO model
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

% Sebastien Besson, Jun 2012


% Input check
ip=inputParser;
ip.addRequired('movieData',@(x) isa(x,'MovieData'));
ip.addRequired('movieInfo',@isstruct);
ip.parse(movieData,movieInfo);

image = movieData.getImage;
updateService = movieData.getSession().getUpdateService();

progressText(0,'Exporting ROIs')
for t=1:size(movieInfo,1)
    
    for i = 1:size(movieInfo(t).xCoord,1)

        roi=omero.model.RoiI();
        roi.setImage(image);
        point = omero.model.PointI;
        point.setCx(omero.rtypes.rdouble(movieInfo(t).xCoord(i,1)));
        point.setCy(omero.rtypes.rdouble(movieInfo(t).yCoord(i,1)));
        point.setTheT(omero.rtypes.rint(t-1));
        point.setTheZ(omero.rtypes.rint(0));
        roi.addShape(point);

        % save
        updateService.saveAndReturnObject(roi);
    end

    progressText(t/size(movieInfo,1))
end
