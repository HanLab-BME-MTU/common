function regexpFilesep=getFilesep(path)
% GETFILESEP returns the file separator associated with a given path as a
% regular expression
%
% 
%   path - A string containing the path which should be replaced
%
% Output: 
%   regexpFilesep - 
% Input:The file separator formatted as a regular expression
%   (to be used in regexp-like function
%
%
% Sebastien Besson, 03/2011
%
   
% Find all separators in the path name
pathSep=unique(regexp(path,'/|\','match'));

%Deal with special cases where pathSep is not exactly '/' or '\'
if isempty(pathSep)
    % Linux path may be:
    % 1 - root path '' without separator
    % 2 - home directory path '~'
    isLinux = @(x) isempty(path) || strcmp(path,'~');
    % Windows path may be:
    % 1 - drive letter with colon e.g. C:
    % 2 -  drive letter without colon e.g. H
    isWindow = @(x) ~isempty(regexp(x,'^[A-Z]$','once')) ||...
        ~isempty(regexp(x,'^[A-Z]:$','once'));
    
    if isLinux(path)
        pathSep='/';
    elseif isWindow(path)
        pathSep='\';
    else
        errordlg(['Error!! Cannot identify the nature of path: ' path]);
        return
    end
elseif numel(pathSep)>1
    errordlg(['Error!! OS conflict in path: ' path]);
    return
else pathSep=pathSep{1};
end

% Return the file separator as a regular expression
regexpFilesep=regexptranslate('escape',pathSep);
end
