function importExternalSegmentation(movieData,varargin)
% importExternalSegmentation imports external masks into the movie infrastructure
%
% This function copies all the folders defined by the InputData field of
% the input parameters under a folder and registers this output in the
% external process.
%
%     importExternalSegmentation(movieData) runs the external segmentation
%     process on the input movie
%
%     importExternalSegmentation(movieData, paramsIn) additionally takes
%     the input parameters
%
%     paramsIn should be a structure with inputs for optional parameters.
%     The parameters should be stored as fields in the structure, with the
%     field names and possible values as described below
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> character string)
%       Optional. A character string specifying the directory to save the
%       masks to. External masks for different channels will be copied
%       under this directory.
%
%       ('InputData' -> Positive integer scalar or vector)
%       Optional. A nChanx1 cell array containing the paths of the folders
%       of the input masks.
%
% Sebastien Besson Nov 2014

importExternalData(movieData, 'ExternalSegmentationProcess', varargin{:});
