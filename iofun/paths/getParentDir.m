% [cpath] = getParentDir(cpath) Returns the parent directory's path from the input path

% Francois Aguet, 11/05/2010

function cpath = getParentDir(cpath)

fsIdx = regexp(cpath, filesep);

if strcmp(cpath(end), filesep)
    cpath = cpath(1:fsIdx(end-1));
else
    cpath = cpath(1:fsIdx(end));
end