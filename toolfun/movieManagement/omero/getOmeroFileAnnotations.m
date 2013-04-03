function [annotations] = getOmeroFileAnnotations(session, imageIDs)


% Input check
ip = inputParser;
ip.addRequired('imageIDs', @isvector);
ip.parse(imageIDs);

namespace = 'hms-tracking';

% Load existing file annotations
% userId = session.getAdminService().getEventContext().userId;
% options = omero.sys.ParametersI;
% options.exp(omero.rtypes.rlong(userId)); %load the annotation for a given user.
metadataService = session.getMetadataService();

% imageIds = java.util.ArrayList;
imageIds = toJavaList(imageIDs, 'java.lang.Long');
annotationTypes = toJavaList('FileAnnotation');
annotatorIds = java.util.ArrayList();

userId = session.getAdminService().getEventContext().userId;
parameters = omero.sys.ParametersI;
parameters.exp(rlong(userId)); %load the annotation for a given user.

annSet = metadataService.loadAnnotations('Image',...
    imageIds, annotationTypes,annotatorIds, parameters);

% Aggregate all annotations into a java.util.ArrayList
annotationList = java.util.ArrayList();
i = annSet.values.iterator;
while (i.hasNext())
    j = i.next().iterator();
    while (j.hasNext())
        annotationList.add(j.next());
    end
end

% Convert java.util.ArrayList into a Matlab array
annotations = toMatlabList(annotationList);