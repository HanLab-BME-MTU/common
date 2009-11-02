% dirList calls 'dir' and returns visible directories.

% Francois Aguet, 11/02/2009
function d = dirList(dpath)

d = dir(dpath);

% remove entries that are not directories
d([d.isdir]==0) = [];

% remove invisible directories
d(cellfun(@isInvisible, {d(:).name})) = [];

function v = isInvisible(directory)
v = strcmp(directory(1), '.');