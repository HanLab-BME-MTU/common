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

if mod(length(varargin),2)~=0
    error('Optional arguments need to be entered as pairs.');
end

if nargin >= 4
    lambda = name2wavelength(lambda);    
    
    % Defaults use values corresponding to optimal imaging conditions
    p.ti0 = 0; % the working distance has no effect under ideal conditions
    p.ni0 = 1.518;
    p.ni = 1.518;
    p.tg0 = 0.17e-3;
    p.tg = 0.17e-3;
    p.ng0 = 1.515;
    p.ng = 1.515;
    p.ns = 1.00;
    p.lambda = lambda;
    p.M = M;
    p.NA = NA;
    p.alpha = asin(p.NA/p.ni);
    p.pixelSize = pixelSize;
else
    p = varargin{1};
end

ru = 8;
psf = vectorialPSF(0,0,0,0,ru,p);

[pG, ~, ~, res] = fitGaussian2D(psf, [0 0 max(psf(:)) 1 0], 'As');
sigma = pG(4);


%===============
% Display
%===============
idx = find(strcmpi(varargin, 'Display'));
if ~isempty(idx) && strcmpi(varargin{idx+1}, 'on')
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
end
