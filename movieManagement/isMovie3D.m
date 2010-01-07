function status = isMovie3D(movieData)

% status = isMovie3D(movieData)
%
% This function returns true if the input movieData structure describes a
% valid 3D movie (as created using setup3DMovieData.m) and false otherwise.
% 
% 
% Hunter Elliott, 11/2009
% 

status = false;

if isfield(movieData,'imSize') && size(movieData.imSize,1) == 3 && ...
        isfield(movieData,'pixelSize_nm') && length(movieData.pixelSize_nm) == 2;
    status = true;
end