function r = bfGetReader(id,varargin)
% Get the reader associated with a given dataset using bioformats
% 
% SYNOPSIS r=bfGetReader(path)
%
% Input 
%    id - the  path of a proprietary file
%
%    debug - Optional. A boolean toggling the logging. Default: true; 
%
% Output
%
%    r - a reader object of class loci.formats.ChannelSeparator
%
% Sebastien Besson, Dec 2011
% Adapted from bfopen.m

% Input check
ip=inputParser;
ip.addRequired('id',@ischar);
ip.addOptional('debug',true,@islogical);
ip.parse(id,varargin{:});
debug=ip.Results.debug;

% load the Bio-Formats library into the MATLAB environment
if ~ismember(which('loci_tools.jar'),javaclasspath)
    javaaddpath(which('loci_tools.jar'));
end

% set LuraWave license code, if available
if exist('lurawaveLicense')
    path = fullfile(fileparts(mfilename('fullpath')), 'lwf_jsdk2.6.jar');
    javaaddpath(path);
    java.lang.System.setProperty('lurawave.license', lurawaveLicense);
end


% initialize logging
if debug
    loci.common.DebugTools.enableLogging('INFO');
else
    loci.common.DebugTools.enableLogging('OFF');
end

r = loci.formats.ChannelFiller();
r = loci.formats.ChannelSeparator(r);

tic
r.setMetadataStore(loci.formats.MetadataTools.createOMEXMLMetadata());
r.setId(id);
