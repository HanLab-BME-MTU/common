function gauss=GaussMask3D(sigma,fSze,cent,cnorm);
% GaussMask3D	create a gaussian 3D mask
%
%    gm=GaussMask3D(sigma,fSze);
%
%    INPUT: sigma  of gauss mask [sigmaX sigmaY sigmaZ]
%           fSze   size of the gauss mask [sizeX sizeY sizeZ]
%                  (odd size required for symmetric mask!)
%           cent   (optional)3D vector with center position [0 0 0] is
%                  center of fSze (default=[0 0 0])
%           cnorm  (optional) select normalization method:
%                  =0 no normalization
%                  =1 (default) norm by total intensity
%
%    OUTPUT: gauss   gaussian mask

% c: 8/05/01 dT

if nargin==2
    cent=[0 0 0];
end;
if nargin<4
    cnorm=1;
end;

gauss=zeros(fSze);
x=([-floor(fSze(1)/2):floor(fSze(1)/2)]-cent(1))./sigma(1);
y=([-floor(fSze(2)/2):floor(fSze(2)/2)]-cent(2))./sigma(2);
z=([-floor(fSze(3)/2):floor(fSze(3)/2)]-cent(3))./sigma(3);
ex=exp(-1/2*(x.^2));
ey=exp(-1/2*(y.^2));
ez=exp(-1/2*(z.^2));
exy=ex'*ey;
exyz=exy(:)*ez;
switch cnorm
case 0
    gauss(:)=exyz(:);    
case 1
    gauss(:)=exyz(:)/sum(exyz(:));
otherwise
    gauss(:)=exyz(:);
end;