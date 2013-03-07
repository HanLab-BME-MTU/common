function [files] = getOmeroFileAnnotations(session, imageIDs)


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
% imageIds.add(java.lang.Long(MD.getImage().getId().getValue))
% annotations = metadataService.loadAnnotation('omero.model.Image', imageIds,...
%     java.util.ArrayList,java.util.ArrayList,omero.sys.ParametersI());

annotationTypes = java.util.ArrayList;
annotationTypes.add('FileAnnotation');

annotatorIds = java.util.ArrayList();
parameters = omero.sys.ParametersI; 
annSet = metadataService.loadAnnotations('Image',...
    toJavaList(imageIDs, 'java.lang.Long'),...
    annotationTypes,annotatorIds, parameters);

itr = annSet.values.iterator;
files = cell(numel(imageIDs), 1);
i = 0;
while (itr.hasNext())
    I = itr.next();
    j = I.iterator();
    i = i + 1;
    while (j.hasNext())
        nextAnn = j.next();
        if strcmp(nextAnn.getNs.getValue, namespace);
            files{i} = nextAnn.getFile;
        end
    end
end