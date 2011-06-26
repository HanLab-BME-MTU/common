% getGaussianPSFsigma returns the standard deviation of the Gaussian approximation of an ideal PSF.
%
% INPUTS     NA        : numerical aperture of the objective
%            M         : magnification of the objective
%            pixelSize : physical pixel size of the CCD in [m]
%            lambda    : emission maximum wavelength of the fluorophore in [m]
%                        -or- fluorophore name
%          {'Display'} : Display PSF and its Gaussian approximation. Optional, default 'off'.
%
% Alternative input: p : parameter structure for PSF calculations (see vectorialPSF.cpp).

% Francois Aguet, October 2010 (Last modified April 6 2011)

function sigma = getGaussianPSFsigma(NA, M, pixelSize, lambda, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('NA', @isscalar);
ip.addRequired('M', @isscalar);
ip.addRequired('pixelSize', @isscalar);
ip.addRequired('lambda', @(x) ischar(x) | isscalar(x))
ip.addParamValue('Display', 'on', @(x) strcmpi(x, 'on') | strcmpi(x, 'off'));
ip.parse(NA, M, pixelSize, lambda, varargin{:});

lambda = name2wavelength(lambda);

% Defaults use values corresponding to optimal imaging conditions
p.ti0 = 0; % working distance has no effect under ideal conditions
p.ni0 = 1.518;
p.ni = 1.518;
p.tg0 = 0.17e-3;
p.tg = 0.17e-3;
p.ng0 = 1.515;
p.ng = 1.515;
p.ns = 1.33;
p.lambda = lambda;
p.M = M;
p.NA = NA;
p.alpha = asin(p.NA/p.ni);
p.pixelSize = pixelSize;


ru = 8;
psf = vectorialPSF(0,0,0,0,ru,p);

[pG, ~, ~, res] = fitGaussian2D(psf, [0 0 max(psf(:)) 1 0], 'As');
sigma = pG(4);


%===============
% Display
%===============
if strcmpi(ip.Results.Display, 'on')
    xa = (-ru+1:ru-1)*p.pixelSize/p.M*1e9;
    
    figure;
    subplot(1,2,1);
    imagesc(xa,xa,psf); colormap(gray(256)); axis image; colorbar;
    title('PSF');
    xlabel('x [nm]');
    ylabel('y [nm]');
    
    subplot(1,2,2);
    imagesc(xa,xa, psf+res.data); colormap(gray(256)); axis image; colorbar;
    title('Gaussian approximation');
    xlabel('x [nm]');
    ylabel('y [nm]');
    linkaxes;
end
