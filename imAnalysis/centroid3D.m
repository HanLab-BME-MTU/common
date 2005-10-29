function ce = centroid3D(img)
% CENTROID compute the centroid of a gray value patch
%
% SYNOPSIS ce = centroid3D(img)
%
% INPUT img : an image 3D patch matrix
% 
% OUTPUT ce : vector with the centroid coordinates (in image coords!)
%
% REMARKS : image patches masked with NaNs will not be counted

% c 19/04/00

s = size(img);
% reshape image only once
img = img(:);
cx=0;
cy=0;
cz=0;
% use centroid so that we get xyz in image coordinates
[x,y,z] = meshgrid(1:s(1),1:s(2),1:s(3));
% nansum in case there are masked regions
cx = nansum(x(:).*img);
cy = nansum(y(:).*img);
cz = nansum(z(:).*img);

ce = [cx, cy, cz]/nansum(img);

% old code (replaced 10/05 by jonas)
%
% for l = 1 : s(2)
%    cx = cx + sum(sum(img(:,l,:).*l));
% end;
% for l = 1 : s(1)
%    cy = cy + sum(sum(img(l,:,:).*l));
% end;
% for l = 1 : s(3)
%    cz = cz + sum(sum(img(:,:,l).*l));
% end;
% sTot=nansum(img(:));
% ce=[cx cy cz]/sTot;