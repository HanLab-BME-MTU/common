function patch=stamp3d(data,patchSize,center)
%STAMP3D copy 3D sub-patch out of larger 3D data set
%
% SYNOPSIS patch=stamp3d(data,patchSize,center)
%
% INPUT data   : 3D data
%       patchSize  : size of patch
%       center : center of patch
%
% OUTPUT patch : 3D patch 

% c: 18/6/01	dT

ds=size(data);
hl=floor(patchSize/2);
hx=min([center(1)-1,hl(1),ds(1)-center(1)]);
hy=min([center(2)-1,hl(2),ds(2)-center(2)]);
hz=min([center(3)-1,hl(3),ds(3)-center(3)]);
patch=data(center(1)-hx:center(1)+hx,center(2)-hy:center(2)+hy,center(3)-hz:center(3)+hz);

