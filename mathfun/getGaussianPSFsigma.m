% getGaussianPSFsigma returns the standard deviation of the Gaussian approximation of an ideal PSF.
%
% INPUTS     NA        : numerical aperture of the objective
%            M         : magnification of the objective
%            pixelSize : physical pixel size of the CCD in [m]
%            lambda    : emission maximum wavelength of the fluorophore in [m]
%                        -or- fluorophore name
%            {display} : Display PSF and its Gaussian approximation. Optional, default 'off'.
%
% Alternative input: p : parameter structure for PSF calculations (see psfVectorial).
%
% For the derivation of the equations used in this function, see Aguet F., PhD thesis, 2009.

% Francois Aguet, October 2010

function sigma = getGaussianPSFsigma(varargin)

if nargin >= 4
    [NA, M, pixelSize, lambda] = deal(varargin{1:4});
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
    p.sf = 3;
    p.mode = 1;
else
    p = varargin{1};
end

ru = 8;
p.k0 = 2*pi/p.lambda;


if nargin==2 || nargin==6
    switch varargin{nargin}
        case 'vectorial'
            integral = psfVectorialRad(0, 0, ru, p);
        case 'scalar'
            integral = psfScalarRad(0, 0, ru, p);
    end
else
    integral = psfVectorialRad(0, 0, ru, p);
end


xyMax = ((2*ru-1)*p.sf-1)/2; % must be fine scale
[x,y] = ndgrid(-xyMax:xyMax);
r = sqrt(x.^2+y.^2);

psf = interp1(0:length(integral)-1, integral, r);
if (p.mode == 1)
    idx = (2*xyMax+1)/p.sf;
    psf = arrayfun(@(k) sum(psf(:,(k-1)*p.sf+1:k*p.sf), 2), 1:idx, 'UniformOutput', false);
    psf = horzcat(psf{:});
    psf = arrayfun(@(k) sum(psf((k-1)*p.sf+1:k*p.sf,:), 1), 1:idx, 'UniformOutput', false);
    psf = vertcat(psf{:});
end

[pG, ~, ~, res] = fitGaussian2D(psf, [0 0 1 1 0], 'As');
sigma = pG(4);

%===============
% Display
%===============
if nargin==5
    xa = (-ru+1:ru-1)*p.pixelSize/p.M*1e9;
    
    figure;
    subplot(1,2,1);
    imagesc(xa,xa,psf); colormap(gray(256)); axis image;
    title('PSF');
    xlabel('x [nm]');
    ylabel('y [nm]');
    
    subplot(1,2,2);
    imagesc(xa,xa, psf+res.data); colormap(gray(256)); axis image;
    title('Gaussian approximation');
    xlabel('x [nm]');
    ylabel('y [nm]');
end




% Bessel functions J0(x) and J1(x)
% Uses the polynomial approximations on p. 369-70 of Abramowitz & Stegun (1972).
% The error in J0 is assumed to be less than or equal to 5 x 10^-8.
function r = J0(x)

% Constants for Bessel function approximation
j0c = [1, -2.2499997, 1.2656208, -0.3163866, 0.0444479, -0.0039444, 0.0002100];
t0c = [-.78539816, -.04166397, -.00003954, 0.00262573, -.00054125, -.00029333, .00013558];
f0c = [.79788456, -0.00000077, -.00552740, -.00009512, 0.00137237, -0.00072805, 0.00014476];

if (x < 0.0)
    x = -x;
end
if (x <= 3.0)
    y = x^2/9.0;
    r = j0c(1) + y*(j0c(2) + y*(j0c(3) + y*(j0c(4) + y*(j0c(5) + y*(j0c(6) + y*j0c(7))))));
else
    y = 3.0/x;
    theta0 = x + t0c(1) + y*(t0c(2) + y*(t0c(3) + y*(t0c(4) + y*(t0c(5) + y*(t0c(6) + y*t0c(7))))));
    f0 = f0c(1) + y*(f0c(2) + y*(f0c(3) + y*(f0c(4) + y*(f0c(5) + y*(f0c(6) + y*f0c(7))))));
    r = sqrt(1.0/x)*f0*cos(theta0);
end



function r = J1(x)

% Constants for Bessel function approximation
j1c = [0.5, -0.56249985, 0.21093573, -0.03954289, 0.00443319, -0.00031761, 0.00001109];
f1c = [0.79788456, 0.00000156, 0.01659667, 0.00017105, -0.00249511, 0.00113653, -0.00020033];
t1c = [-2.35619449, 0.12499612, 0.00005650, -0.00637897, 0.00074348, 0.00079824, -0.00029166];

sign = 1.0;
if (x < 0.0)
    x = -x;
    sign = -sign;
end
if (x <= 3.0)
    y = x^2/9.0;
    r = x*(j1c(1) + y*(j1c(2) + y*(j1c(3) + y*(j1c(4) + y*(j1c(5) + y*(j1c(6) + y*j1c(7)))))));
else
    y = 3.0/x;
    theta1 = x + t1c(1) + y*(t1c(2) + y*(t1c(3) + y*(t1c(4) + y*(t1c(5) + y*(t1c(6) + y*t1c(7))))));
    f1 = f1c(1) + y*(f1c(2) + y*(f1c(3) + y*(f1c(4) + y*(f1c(5) + y*(f1c(6) + y*f1c(7))))));
    r = sqrt(1.0/x)*f1*cos(theta1);
end
r = sign*r;



function [L dL] = L_theta(theta, p)

ni2sin2theta = p.ni^2*sin(theta)^2;
groot = sqrt(p.ng^2 - ni2sin2theta);
g0root = sqrt(p.ng0^2 - ni2sin2theta);
i0root = sqrt(p.ni0^2 - ni2sin2theta);

ci = p.ni*(p.tg0/p.ng0 + p.ti0/p.ni0 - p.tg/p.ng);

L = p.ni*ci*cos(theta) + p.tg*groot - p.tg0*g0root - p.ti0*i0root;
dL = p.ni*sin(theta) * (-ci + p.ni*cos(theta)*(p.tg0/g0root + p.ti0/i0root - p.tg/groot));



function integral = psfScalarRad(xp, yp, ru, p)

xystep = p.pixelSize/p.M;
NA2 = p.NA*p.NA;
Ni02 = p.ni0*p.ni0;
Ni2 = p.ni*p.ni;
Ng02 = p.ng0*p.ng0;
Ng2 = p.ng*p.ng;
A0 = Ni2*Ni2/(NA2*NA2);

% constant component of OPD
ci = p.ni*(p.tg0/p.ng0 + p.ti0/p.ni0 - p.tg/p.ng);

xp_n = xp*p.sf/xystep;
yp_n = yp*p.sf/xystep;
rn = 1 + floor(sqrt(xp_n*xp_n + yp_n*yp_n));

xyMax = ((2*ru-1)*p.sf-1)/2; % must be fine scale
rMax = ceil(sqrt(2.0)*xyMax) + rn + 1; % +1 for interpolation, dx, dy


integral = zeros(1, rMax);
ud = 3.0*p.sf;

theta0 = p.alpha;
[~, dL] = L_theta(theta0, p);

w_exp = abs(dL); % missing p.k0 !

cst = 0.975;
while (cst >= 0.9)
    [~, dL] = L_theta(cst*theta0, p);
    if (abs(dL) > w_exp)
        w_exp = abs(dL);
    end
    cst = cst - 0.025;
end
w_exp = w_exp*p.k0;

for ri = 0:rMax-1
    
    r = xystep/p.sf*ri;
    constJ = p.k0*r*p.ni; % = w_J;
    
    if (w_exp > constJ)
        nSamples = 4 * floor(1.0 + p.alpha*w_exp/pi);
    else
        nSamples = 4 * floor(1.0 + p.alpha*constJ/pi);
    end
    if (nSamples < 20)
        nSamples = 20;
    end
    step =  p.alpha/nSamples;
    iconst = (step/ud)^2;
    
    % Simpson's rule
    sum_I0 = 0.0;
    
    for n = 1:nSamples/2-1
        theta = 2.0*n*step;
        sintheta = sin(theta);
        costheta = cos(theta);
        ni2sin2theta = p.ni*p.ni*sintheta*sintheta;
        bessel_0 = 2.0*J0(constJ*sintheta)*sintheta*costheta; % 2.0 factor : Simpson's rule
        expW = exp(1i*p.k0*(ci*p.ni*costheta + p.tg*sqrt(Ng2-ni2sin2theta) - p.tg0*sqrt(Ng02-ni2sin2theta) - p.ti0*sqrt(Ni02-ni2sin2theta)));
        sum_I0 = sum_I0 + expW*bessel_0;
    end
    for n = 1:nSamples/2
        theta = (2.0*n-1.0)*step;
        sintheta = sin(theta);
        costheta = cos(theta);
        ni2sin2theta = p.ni*p.ni*sintheta*sintheta;
        bessel_0 = 4.0*J0(constJ*sintheta)*sintheta*costheta;
        expW = exp(1i*p.k0*(ci*p.ni*costheta + p.tg*sqrt(Ng2-ni2sin2theta) - p.tg0*sqrt(Ng02-ni2sin2theta) - p.ti0*sqrt(Ni02-ni2sin2theta)));
        sum_I0 = sum_I0 + expW*bessel_0;
    end
    % theta = alpha;
    bessel_0 = J0(p.k0*r*p.NA)*cos(p.alpha)*sin(p.alpha);
    expW = exp(1i*p.k0*(ci*sqrt(Ni2-NA2) + p.tg*sqrt(Ng2-NA2) - p.tg0*sqrt(Ng02-NA2) - p.ti0*sqrt(Ni02-NA2)));
    sum_I0 = sum_I0 + expW*bessel_0;
    
    integral(ri+1) = A0 * abs(sum_I0)*abs(sum_I0) * iconst;
end



function integral = psfVectorialRad(xp, yp, ru, p)

xystep = p.pixelSize/p.M;

% constant component of OPD
ci = p.ni*(p.tg0/p.ng0 + p.ti0/p.ni0 - p.tg/p.ng);

xp_n = xp*p.sf/xystep;
yp_n = yp*p.sf/xystep;
rn = 1 + floor(sqrt(xp_n*xp_n + yp_n*yp_n));

xyMax = ((2*ru-1)*p.sf-1)/2; % must be fine scale
rMax = ceil(sqrt(2.0)*xyMax) + rn + 1; % +1 for interpolation, dx, dy

integral = zeros(1, rMax);
ud = 3.0*p.sf;

theta0 = p.alpha;
[~, dL] = L_theta(theta0, p);
w_exp = abs(dL); % missing p.k0 !

cst = 0.975;
while (cst >= 0.9)
    [~, dL] = L_theta(cst*theta0, p);
    if (abs(dL) > w_exp)
        w_exp = abs(dL);
    end
    cst = cst - 0.025;
end
w_exp = w_exp * p.k0;

for ri = 0:rMax-1
    
    r = xystep/p.sf*ri;
    constJ = p.k0*r*p.ni; % = w_J;
    
    if (w_exp > constJ)
        nSamples = 4 * floor(1.0 + p.alpha*w_exp/pi);
    else
        nSamples = 4 * floor(1.0 + p.alpha*constJ/pi);
    end
    if (nSamples < 20)
        nSamples = 20;
    end
    
    step =  p.alpha/nSamples;
    iconst = step/ud;
    
    % Simpson's rule
    sum_I0 = 0.0;
    sum_I1 = 0.0;
    sum_I2 = 0.0;
    
    for n = 1:nSamples/2-1
        theta = 2.0*n*step;
        sintheta = sin(theta);
        costheta = cos(theta);
        sqrtcostheta = sqrt(costheta);
        ni2sin2theta = p.ni^2*sintheta*sintheta;
        nsroot = sqrt(p.ns^2 - ni2sin2theta);
        ngroot = sqrt(p.ng^2 - ni2sin2theta);
        
        ts1ts2 = 4.0*p.ni*costheta*ngroot;
        tp1tp2 = ts1ts2;
        tp1tp2 = tp1tp2 / ((p.ng*costheta + p.ni/p.ng*ngroot) * (p.ns/p.ng*ngroot + p.ng/p.ns*nsroot));
        ts1ts2 = ts1ts2 / ((p.ni*costheta + ngroot) * (ngroot + nsroot));
        
        bessel_0 = 2.0*J0(constJ*sintheta)*sintheta*sqrtcostheta; % 2.0 factor : Simpson's rule
        bessel_1 = 2.0*J1(constJ*sintheta)*sintheta*sqrtcostheta;
        if (constJ ~= 0.0)
            bessel_2 = 2.0*bessel_1/(constJ*sintheta) - bessel_0;
        else
            bessel_2 = 0.0;
        end
        bessel_0 = bessel_0 * (ts1ts2 + tp1tp2/p.ns*nsroot);
        bessel_1 = bessel_1 * (tp1tp2*p.ni/p.ns*sintheta);
        bessel_2 = bessel_2 * (ts1ts2 - tp1tp2/p.ns*nsroot);
        
        expW = exp(1i*p.k0*(ci*p.ni*costheta + p.tg*ngroot - p.tg0*sqrt(p.ng0^2-ni2sin2theta) - p.ti0*sqrt(p.ni0^2-ni2sin2theta)));
        sum_I0 = sum_I0 + expW*bessel_0;
        sum_I1 = sum_I1 + expW*bessel_1;
        sum_I2 = sum_I2 + expW*bessel_2;
    end
    for n = 1:nSamples/2
        theta = (2.0*n-1)*step;
        sintheta = sin(theta);
        costheta = cos(theta);
        sqrtcostheta = sqrt(costheta);
        ni2sin2theta = p.ni^2*sintheta*sintheta;
        nsroot = sqrt(p.ns^2 - ni2sin2theta);
        ngroot = sqrt(p.ng^2 - ni2sin2theta);
        
        ts1ts2 = 4.0*p.ni*costheta*ngroot;
        tp1tp2 = ts1ts2;
        tp1tp2 = tp1tp2 / ((p.ng*costheta + p.ni/p.ng*ngroot) * (p.ns/p.ng*ngroot + p.ng/p.ns*nsroot));
        ts1ts2 = ts1ts2 / ((p.ni*costheta + ngroot) * (ngroot + nsroot));
        
        bessel_0 = 4.0*J0(constJ*sintheta)*sintheta*sqrtcostheta; % 4.0 factor : Simpson's rule
        bessel_1 = 4.0*J1(constJ*sintheta)*sintheta*sqrtcostheta;
        if (constJ ~= 0.0)
            bessel_2 = 2.0*bessel_1/(constJ*sintheta) - bessel_0;
        else
            bessel_2 = 0.0;
        end
        bessel_0 = bessel_0 * (ts1ts2 + tp1tp2/p.ns*nsroot);
        bessel_1 = bessel_1 * (tp1tp2*p.ni/p.ns*sintheta);
        bessel_2 = bessel_2 * (ts1ts2 - tp1tp2/p.ns*nsroot);
        
        expW = exp(1i*p.k0*(ci*p.ni*costheta + p.tg*ngroot - p.tg0*sqrt(p.ng0^2-ni2sin2theta) - p.ti0*sqrt(p.ni0^2-ni2sin2theta)));
        sum_I0 = sum_I0 + expW*bessel_0;
        sum_I1 = sum_I1 + expW*bessel_1;
        sum_I2 = sum_I2 + expW*bessel_2;
    end
    % theta = alpha;
    sintheta = sin(p.alpha);
    costheta = cos(p.alpha);
    sqrtcostheta = sqrt(costheta);
    nsroot = sqrt(p.ns^2 - p.NA^2);
    ngroot = sqrt(p.ng^2 - p.NA^2);
    
    ts1ts2 = 4.0*p.ni*costheta*ngroot;
    tp1tp2 = ts1ts2;
    tp1tp2 = tp1tp2 / ((p.ng*costheta + p.ni/p.ng*ngroot) * (p.ns/p.ng*ngroot + p.ng/p.ns*nsroot));
    ts1ts2 = ts1ts2 / ((p.ni*costheta + ngroot) * (ngroot + nsroot));
    
    bessel_0 = J0(constJ*sintheta)*sintheta*sqrtcostheta;
    bessel_1 = J1(constJ*sintheta)*sintheta*sqrtcostheta;
    if (constJ ~= 0.0)
        bessel_2 = 2.0*bessel_1/(constJ*sintheta) - bessel_0;
    else
        bessel_2 = 0.0;
    end
    bessel_0 = bessel_0 * (ts1ts2 + tp1tp2/p.ns*nsroot);
    bessel_1 = bessel_1 * (tp1tp2*p.ni/p.ns*sintheta);
    bessel_2 = bessel_2 * (ts1ts2 - tp1tp2/p.ns*nsroot);
    
    expW = exp(1i*p.k0*(ci*sqrt(p.ni^2-p.NA^2) + p.tg*ngroot - p.tg0*sqrt(p.ng0^2-p.NA^2) - p.ti0*sqrt(p.ni0^2-p.NA^2)));
    sum_I0 = sum_I0 + expW*bessel_0;
    sum_I1 = sum_I1 + expW*bessel_1;
    sum_I2 = sum_I2 + expW*bessel_2;
    
    sum_I0 = abs(sum_I0);
    sum_I1 = abs(sum_I1);
    sum_I2 = abs(sum_I2);
    
    integral(ri+1) = 8.0*pi/3.0*real(sum_I0*sum_I0 + 2.0*sum_I1*sum_I1 + sum_I2*sum_I2) * iconst*iconst;
end