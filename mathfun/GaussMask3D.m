function gauss=GaussMask3D(sigma,fSze,cent,cnorm)
% GaussMask3D	create a gaussian 3D mask
%
%    SZNOPSIS gauss =GaussMask3D(sigma,fSze,cent,cnorm);
%
%    INPUT: sigma  of gauss mask [sigmaX sigmaY sigmaZ]
%           fSze   size of the gauss mask [sizeX sizeY sizeZ]
%                  (odd size required for symmetric mask!)
%           cent   (optional)3D vector with center position [0 0 0] is
%                  center of fSze (default=[0 0 0])
%           cnorm  (optional) select normalization method:
%                  =0 (default) no normalization - max of gauss will be 1
%                  =1 norm so that integral will be 1 (not quite, as the
%                  gaussMask can never be infinitely large
%                      
%
%    OUTPUT: gauss   gaussian mask

% c: 8/05/01 dT

% test input
if nargin==2
    cent=[0 0 0];
end;
if nargin<4
    cnorm=0;
end;
if ~all((fSze-1)/2 == floor((fSze-1)/2))
    error('gauss mask needs odd size!')
end

% to get the accurate pixel intensities, use the cumsum of the gaussian
% (taken from normcdf.m). As the origin is in the center of the center
% pixel, the intensity of pixel +2 is the integral from 1.5:2.5, which,
% fortunately, has the spacing 1.

gauss=zeros(fSze);
x=([-fSze(1)/2:fSze(1)/2]-cent(1))./sigma(1);
y=([-fSze(2)/2:fSze(2)/2]-cent(2))./sigma(2);
z=([-fSze(3)/2:fSze(3)/2]-cent(3))./sigma(3);
ex = diff(0.5 * erfc(-x./sqrt(2)));
ey = diff(0.5 * erfc(-y./sqrt(2)));
ez = diff(0.5 * erfc(-z./sqrt(2)));

% construct the 3D matrix (nice work by Dom!)
exy=ex'*ey;
exyz=exy(:)*ez;

% norm Gauss
switch cnorm
case 0 % maximum of Gauss has to be 1
    gauss(:) = exyz(:)*((2*pi)^1.5*sigma(1)*sigma(2)*sigma(3));  
    
case 1 % the whole erfc thing is already normed, so nothing to do here.  
    gauss(:) = exyz(:);
    
otherwise
    gauss(:) = exyz(:);
    
end