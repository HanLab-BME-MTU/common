function data=pushstamp3d(data,patch,center)
%PUSHSTAMP3D writes a 3Dsubimage into a larger image 
%
% SYNOPSIS data=pushstamp3d(data,patchSize,center)
%
% INPUT data   : 3D data
%            patchSize  : size of patch
%            center : center of patch (in 3D data)
%
% OUTPUT data : 3D data 

% c: 18/6/01	dT

ds=size(data);
hl=floor(size(patch)/2);
hx=min([center(1)-1,hl(1),ds(1)-center(1)]);
hy=min([center(2)-1,hl(2),ds(2)-center(2)]);
hz=min([center(3)-1,hl(3),ds(3)-center(3)]);
data(center(1)-hx:center(1)+hx,center(2)-hy:center(2)+hy,center(3)-hz:center(3)+hz)=patch(hl(1)+1-hx:hl(1)+1+hx,hl(2)+1-hy:hl(2)+1+hy,hl(3)+1-hz:hl(3)+1+hz);

