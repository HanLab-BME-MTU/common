% Francois Aguet, November 2010

function dirName = getDirFromPath(dpath)

idx = regexp(dpath, filesep);
if idx(end) == length(dpath)
    dirName = dpath(idx(end-1)+1:end-1);
else
    dirName = dpath(idx(end)+1:end);
end