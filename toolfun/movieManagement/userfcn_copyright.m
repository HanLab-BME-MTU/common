function info = userfcn_copyright (type)

% This function returns the information of LCCB software
%
% Input: 
%
%   type - The type of information returned
%       ('year' - release year)
%       ('version' - version number)
%       ('all' - output is in the format e.g. '© 2011 LCCB   Version 3.0')
%
% Chuangang Ren
% 11/2010

str_year = '2011';
str_version = '2.0';

if nargin < 1
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


