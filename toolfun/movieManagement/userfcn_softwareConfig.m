function [info openHelpFile] = userfcn_softwareConfig(handles, type)
%
% This is a user-defined function used in LCCB software. 
% It is called when any GUI is generated. It configures the software
% information.
% Now it is able to set the following parameters for the software
% 1. copyright and software version information displayed at the upper left
%    of each GUI
% 2. Whether to open original help file when user clicks on Help icon in
%    the GUI
%
% Input: 
%
%   handles - reserved, to be used for distinguishing different GUIs
%             in a future version
%   type - The type of information returned
%       ('year' - release year)
%       ('version' - version number)
%       ('all' - output is in the format e.g. 'Copyright 2011 LCCB   Version 3.0')
%
% Output:
%
%   info - String: copyright and version information
%   openHelpFile - 1 or 0, whether to open original help file when user
%                  clicks on Help icon in the GUI. 1 - yes, 0 - no.
%
% Chuangang Ren
% 11/2010

% Set 1 or 0 
openHelpFile = 1;

% Set year and version information
str_year = '2011';
str_version = '2.1';

if nargin < 2
    type = 'all';
end

if ~ischar(type)
   error('User-defined: Input is not a char.') 
end

switch lower(type)
    
    case 'year'
        info = str_year;
        
    case 'version'
        info = str_version;
        
    case 'all'
        info = sprintf('Copyright %s LCCB    Version %s', str_year, str_version);
        
    otherwise
        
end

