function [data] = locateSegmentation(data);
% locate segmentation data in the specified paths, and store the location 
% as fields in the data structure
%
% INPUT     data    :    experiment structure, which has to contain a field
%                       .source
%                        which is the path to the data location
% OUTPUT    two new fields are added to (or overwritten in) data, which are
%                   data.segmentDataFilePath
%                   data.segmentDataFileName
%
% NOTE: If the folder containing the segmentation data has a standardized
% name (e.g. if it's always called 'SegmentImages' or something along these
% lines), the function can be adapted to be more user-friendly by
% automatically setting the path to the standard-name subfolder, instead of
% just the source path - this would mean less user clicking
%
% Dinah Loerke, last modified April 20, 2008

od = cd;

% loop over all entries in the structure to enter the image data necessary
% for the detection input
for i=1:length(data)
    
    path = data(i).source;
    % change to source directory
    cd(path);
    
    [oriImageName, oriImagePath] = uigetfile('.tif',['Select first segmentation image']); 
    
    data(i).segmentDataFileName = oriImageName;
    data(i).segmentDataFilePath = oriImagePath;
        
    cd(od);
end


end % of function