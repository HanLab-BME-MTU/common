function ce = centroid3D(img)
% CENTROID compute the centoid of a gray value patch
%
%
% SYNOPSIS ce = centroid3D(img)
%
% INPUT img : an image 3D patch matrix
% 
% OUTPUT ce : vector with the centroid coordinates

% c 19/04/00

s = size(img);
cx=0;
cy=0;
cz=0;
for l = 1 : s(2)
   cx = cx + sum(sum(img(:,l,:).*l));
end;
for l = 1 : s(1)
   cy = cy + sum(sum(img(l,:,:).*l));
end;
for l = 1 : s(3)
   cz = cz + sum(sum(img(:,:,l).*l));
end;
sTot=sum(img(:));
ce=[cx cy cz]/sTot;