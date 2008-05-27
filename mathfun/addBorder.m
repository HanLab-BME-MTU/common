function array = addBorder(array,border)
%ADDBORDER adds a mirrored border around 1-3d arrays
%
% SYNOPSIS: array = addBorder(array,border)
%
% INPUT array: d-dimensional array, where d can be 1-3
%		border : 1-by-d array indicating by how many pixels the array
%                should be extended. Usually, that is going to be half of
%                your array size.
%                Pass a 2-by-d array if there is unequal padding on the two
%                sides of the array
%
% OUTPUT array: array with added borders
%
% REMARKS addBorder reflects the array to pad it
%
% created with MATLAB ver.: 7.6.0.324 (R2008a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 26-May-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
if nargin < 2 || isempty(array) || isempty(border)
    error('please supply at least two non-empty input arguments')
end

nDims = ndims(array);
if nDims == 1
    array = array(:);
end
if nDims > 3
    error('please supply maximum 3D arrays')
end
if size(border,2) ~= nDims
    error('border needs a column for every dimension')
end
if size(border,1) == 1
    border = [border;border];
end

% flip along first dimension
array = cat(1,flipdim(array(1:border(1,1),:,:),1),...
    array,flipdim(array(end-border(1,2)+1:end,:,:),1));
if nDims > 1
    % flip along second dimension
    array = cat(2,flipdim(array(:,1:border(1,2),:),2),...
    array,flipdim(array(:,end-border(2,2)+1:end,:),2));
end

if nDims > 2
    % flip along third dimension
     array = cat(3,flipdim(array(:,:,1:border(1,3)),3),...
    array,flipdim(array(:,:,end-border(2,3)+1:end),3));
end