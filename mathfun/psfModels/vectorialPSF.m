% h = vectorialPSF(xp, yp, zp, z, ru, p)
%
% Computes microscope PSFs using the model described in:
% Aguet et al., Opt. Express 17(8), 6829-6848 (2009)
%
%   INPUTS:
%   xp    : x-position of the source in object space coordinates.
%   yp    : y-position of the source.
%   zp    : z-position of the source.
%   z     : vector of focus/acquisition positions.
%   ru    : radial size of the psf, in pixels. The resulting psf has dimensions [2*ru-1 2*ru-1 length(z)].
%
%   p     : parameter structure for the system properties, with the fields
%            Ti0       : working distance of the objective
%            Ni0       : immersion medium refractive index, design value
%            Ni        : immersion medium refractive index, experimental value
%            Tg0       : coverslip thickness, design value
%            Tg        : coverslip thickness, experimental value
%            Ng0       : coverslip refractive index, design value
%            Ng        : coverslip refractive index, experimental value
%            Ns        : sample refractive index
%            lambda    : emission wavelength
%            M         : magnification
%            NA        : numerical aperture
%            alpha     : angular aperture, alpha = atan(NA/Ni)
%            pixelSize : physical size (width) of the CCD pixels
%            f         : (optional, default: 3) oversampling factor to approximate CCD integration
%            mode      : (optional, default: 1) if mode = 0, returns oversampled PSF
%
%   All spatial units in [m].
% 
%
%   Example structure for 'p':
%
%           p.Ti0 = 1.9000e-04
%           p.Ni0 = 1.5180
%            p.Ni = 1.5180
%           p.Tg0 = 1.7000e-04
%            p.Tg = 1.7000e-04
%           p.Ng0 = 1.5150
%            p.Ng = 1.5150
%            p.Ns = 1
%        p.lambda = 5.5000e-07
%             p.M = 100
%            p.NA = 1.4500
%         p.alpha = asin(p.NA/p.Ni);
%     p.pixelSize = 6.4500e-06

% Francois Aguet, last modified: 11/01/2010
function vectorialPSF