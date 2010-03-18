function [dFdXc dFdYc dFdA dFds dFdl dFdt] = ...
    dLSegment2DJacobian(xRange, yRange, xC, yC, A, sigmaPSF, l, theta)
% Partial derivatives of a 2-dimensional diffraction-limited segment.
% [dFdXc dFdXy dFdA dFds dFdl dFdt] =
%           dLSegment2D_dFdA(xRange, yRange, xC, xC, A, sigmaPSF, l, theta)
%
% parameters:
% (xRange, yRange)   2 vectors representing the 2-dimensional domain (e.g.
%                    xRange = -10:.1:10, yRange = -5:.1:5
%
% (xC,yC)            center of the segment
%
% A                  amplitude of the segment
%
% sigmaPSF           half width of the gaussian PSF model.
%
% l                  length of the segment
%
% theta              orientation of the segment
%
% output:
% dFdXc dFdXy dFdA dFds dFdl dFdt are NxM matrices where N = numel(xRange)
% and M = numel(yRange). 
%
% Sylvain Berlemont, 2010

ct = cos(theta);
st = sin(theta);

l = l / 2;

[X Y] = meshgrid(xRange, yRange);

X = X - xC;
Y = Y - yC;

dFdXc = (1/2).*A.*exp(1).^((-1/2).*sigmaPSF.^(-2).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).^2).*sigmaPSF.^(-2).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-1).*(exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*(y+(-1).*yc).*sin(theta)).^2).*(2.*pi.^(-1)).^(1/2).*sigmaPSF.*cos(theta)+(-1).*exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)).^2).*(2.*pi.^(-1)).^(1/2).*sigmaPSF.*cos(theta)+(Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)))+Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta)))).*sin(theta).*(((-1).*y+yc).*cos(theta)+(x+(-1).*xc).*sin(theta)));

dFdYc = (1/2).*A.*exp(1).^((-1/2).*sigmaPSF.^(-2).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).^2).*sigmaPSF.^(-2).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-1).*(exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*(y+(-1).*yc).*sin(theta)).^2).*(2.*pi.^(-1)).^(1/2).*sigmaPSF.*sin(theta)+(-1).*exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)).^2).*(2.*pi.^(-1)).^(1/2).*sigmaPSF.*sin(theta)+cos(theta).*(Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)))+Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta)))).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)));

dFdA = (1/2).*exp(1).^((-1/2).*sigmaPSF.^(-2).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).^2).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-1).*(Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)))+Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta))));

dFds = (1/2).*A.*exp(1).^((-1/2).*sigmaPSF.^(-2).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).^2).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-2).*(exp(1).^((-1/2).*l.^2.*sigmaPSF.^(-2)).*l.*(2.*pi.^(-1)).^(1/2).*sigmaPSF.^(-2).*(Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)))+Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta))))+sigmaPSF.^(-3).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).*(Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)))+Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta)))).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).^2+(2.*pi).^(-1/2).*sigmaPSF.^(-2).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).*(exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*(y+(-1).*yc).*sin(theta)).^2).*((-1).*l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta))+(-1).*exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)).^2).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta))));

dFdl = (1/2).*A.*exp(1).^((-1/2).*sigmaPSF.^(-2).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).^2).*(2.*pi).^(-1/2).*sigmaPSF.^(-1).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-2).*((exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)).^2)+exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+(-2).*x.*cos(theta)+2.*xc.*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta)).^2)).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1))+(-2).*exp(1).^((-1/2).*l.^2.*sigmaPSF.^(-2)).*(Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)))+Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+(-2).*x.*cos(theta)+2.*xc.*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta)))));

dFdt = (1/2).*A.*exp(1).^((-1/2).*sigmaPSF.^(-2).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).^2).*Erf(2.^(-1/2).*l.*sigmaPSF.^(-1)).^(-1).*(exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*(y+(-1).*yc).*sin(theta)).^2).*(2.*pi).^(-1/2).*sigmaPSF.^(-1).*((-2).*(y+(-1).*yc).*cos(theta)+2.*(x+(-1).*xc).*sin(theta))+exp(1).^((-1/8).*sigmaPSF.^(-2).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)).^2).*(2.*pi.^(-1)).^(1/2).*sigmaPSF.^(-1).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta))+sigmaPSF.^(-2).*(Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+2.*(x+(-1).*xc).*cos(theta)+2.*(y+(-1).*yc).*sin(theta)))+Erf((1/2).*2.^(-1/2).*sigmaPSF.^(-1).*(l+(-2).*(x+(-1).*xc).*cos(theta)+(-2).*y.*sin(theta)+2.*yc.*sin(theta)))).*((y+(-1).*yc).*cos(theta)+((-1).*x+xc).*sin(theta)).*((x+(-1).*xc).*cos(theta)+(y+(-1).*yc).*sin(theta)));

