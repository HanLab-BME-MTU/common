function out = vectorialPSF(xp, zvec, nx, params, varargin)
%VECTORIALPSF Vectorial (Richards–Wolf) PSF in stratified media (all MATLAB).
%
% out = vectorialPSF(xp, zvec, nx, params)
% out = vectorialPSF(..., 'derivXP', true)  % also returns dPSF/dxp
%
% Inputs
%   xp        : emitter lateral offset [m] (scalar; shift along +x)
%   zvec      : axial positions [m], 1×nz (z = 0 can be coverslip/objective
%               focal reference; sign convention: positive toward camera)
%   nx        : lateral output size (image is nx × nx for each z)
%   params    : struct with fields (SI units)
%        .lambda   emission wavelength in medium above objective [m]
%        .ni       refractive index immersion (objective side)
%        .ns       refractive index sample medium
%        .ng       refractive index coverslip (glass)
%        .NA       numerical aperture
%        .M        magnification (for pixel sampling on camera)
%        .pixel    camera pixel pitch [m]
%        .t_ig     thickness (immersion->glass) [m]
%        .t_gs     thickness (glass->sample) [m]
%        .Ntheta   number of polar (theta) samples (even; default 400)
%        .Nr       number of radial r samples in pupil (uses GL nodes;
%                   default chosen from Ntheta)
%        .sf       oversampling factor within one camera pixel (default 3)
%   Optional name/value:
%        'derivXP' : logical, also compute dPSF/dxp (default false)
%
% Output (struct)
%   .PSF     : nx×nx×nz double
%   .dxp     : nx×nx×nz (only if 'derivXP',true)
%   .meta    : copy of params plus computed constants
%
% Notes
%   - Direct translation of vectorialPSF.cpp via ChatGPT (by Sangyoon)
%   - Uses Richards–Wolf field integrals with Bessel kernels J0,J1,J2.
%   - Stratified s/p transmission at immersion→glass→sample interfaces.
%   - Uses Gauss–Legendre in theta. Vectorized over (x,y,z) grids.
%   - If you only need 2-D PSF at a single z, pass zvec as scalar.
%
%   This is written for clarity & robustness (no GPU, no MEX).
%   If you later want speed, we can fuse loops and add GPU arrays.
%

% ------------------- sanity and defaults -------------------
if ~isfield(params,'Ntheta'), params.Ntheta = 400; end
if mod(params.Ntheta,2)~=0, params.Ntheta = params.Ntheta+1; end
if ~isfield(params,'sf'),     params.sf     = 3;   end
if ~isfield(params,'M'),      params.M      = 100; end

wantDXP = any(strcmpi(varargin,'derivXP')) && any([varargin{2:end}]==true);

lambda = params.lambda;
k0     = 2*pi/lambda;            % vacuum k0 if lambda outside? assume already in medium
ni     = params.ni;
ns     = params.ns;
ng     = params.ng;
NA     = params.NA;

alpha  = asin(min(NA/ni,1));     % maximal pupil angle in immersion
ni2    = ni^2; ns2 = ns^2; ng2 = ng^2;

% camera/grid sampling
dx_cam = params.pixel/params.M/params.sf;   % sampling in object space
x = ((-(nx*params.sf-1)/2):((nx*params.sf-1)/2))*dx_cam;
[X,Y] = meshgrid(x,x);   % object-plane coords
Rho   = hypot(X-xp, Y);  % lateral distance from emitter offset
Phi   = atan2(Y, X-xp);
nz = numel(zvec);

% ------------------- quadrature in theta -------------------
N = params.Ntheta;
% Gauss-Legendre nodes/weights on [0,alpha]
[t, w] = gaussLegendre(0, alpha, N);

ct = cos(t);  st = sin(t);  sct = sqrt(ct);

% Stratified s/p transmission (immersion -> glass -> sample)
% Fresnel at each interface for s and p:
%   t_s = 2 n1 cos(theta1) / (n1 cos(theta1)+n2 cos(theta2))
%   t_p = 2 n1 cos(theta1) / (n2 cos(theta1)+n1 cos(theta2))
% with Snell: n1*sin(theta1) = n2*sin(theta2)
%
% We form total transmission T = t12 * t23 with internal angles.
% Handle TIR by complex angles (use sqrt with complex).
ni2st2 = ni2*(st.^2);
c2g = sqrt(complex(1 - ni2st2/ng2));     % cos(theta_g) up to constant
c2s = sqrt(complex(1 - ni2st2/ns2));     % cos(theta_s)

% interface 1: immersion(ni) -> glass(ng)
ts_ig = 2*ni*ct ./ (ni*ct + ng*c2g);
tp_ig = 2*ni*ct ./ (ng*ct + ni*c2g);

% interface 2: glass(ng) -> sample(ns), angle inside glass is arcsin(ni/ng sin t)
% cos(theta at glass side)
c_glass = c2g;  % already cos(theta_g)
% sin(theta_g) not needed explicitly; use same ni2st2 trick
% For passage glass->sample, need cos in glass and cos in sample:
ts_gs = 2*ng*c_glass ./ (ng*c_glass + ns*c2s);
tp_gs = 2*ng*c_glass ./ (ns*c_glass + ng*c2s);

Ts = ts_ig .* ts_gs;
Tp = tp_ig .* tp_gs;

% Phase through layers (geometric OPD)
% Add thickness if provided (defaults 0)
tig = getfielddef(params,'t_ig',0);
tgs = getfielddef(params,'t_gs',0);

k_ni = k0*ni; k_ng = k0*ng; k_ns = k0*ns;

% axial phase terms for each theta (z will come later)
phi_ig = k_ng*tig.*c_glass + k_ns*tgs.*c2s - k_ni*(tig+tgs).*ct;  % relative vs immersion
phase_static = exp(1i*phi_ig);    % size 1×N

% Polarization weight factors for RW integrals (following Born & Wolf conv.)
% The (a0,a1,a2) correspond to components producing J0,J1,J2 kernels.
% Here we follow Sheppard/Török / Gibson–Lanni style aggregation:
A0 = sct .* ct .* Ts;                    % s-like contribution
A2 = sct .* (1-ct) .* Tp;                % p-like transverse term
A1 = sct .* st .* Tp;                    % p-like longitudinal coupling

% Precompute common pieces across radii: argument to Bessel is k*R*sin(theta)
% We will evaluate per pixel with vectorized BesselJ at those arguments.

% -------------- accumulate fields for all (x,y,z) ----------------
PSF  = zeros(nx*params.sf, nx*params.sf, nz, 'double');
dPSFdxp = [];

% axial grid phase factor exp(i k ns z cos(theta_s)) relative to immersion focus
% We propagate in sample with ns and in immersion with ni to keep mismatch;
% the commonly used model uses the sample space for emitter—adjust if needed.
for iz = 1:nz
    z = zvec(iz);
    % axial phase relative to focus:
    % choose sample side propagation for physical emitter in sample:
    axial = exp(1i * (k_ns*z).*c2s);
    Wz = phase_static .* axial;      % 1×N
    % kernels at each pixel radius Rho:
    rho = Rho;  % (ny×nx)
    arg = (k_ni * st).';             % column N×1
    % For vectorization, reshape:
    rr  = rho(:).';                  % 1×P
    kr  = arg * rr;                  % (N×P)
    % J0,J1,J2:
    J0 = besselj(0, kr);
    J1 = besselj(1, kr);
    J2 = besselj(2, kr);

    % Phase term e^{i k ni xp cos(phi) sin theta} is already included via shifting Rho by xp
    % so nothing else here.

    % Build integrals (Gauss–Legendre weights)
    wcol = (w(:) .* Wz(:)).';     % 1×N
    I0 = ( (A0(:).'.*wcol) * J0 );        % 1×P
    I1 = ( (A1(:).'.*wcol) * J1 ) .* cos(Phi(:).');  % coupling with cos(phi)
    I2 = ( (A2(:).'.*wcol) * J2 ) .* cos(2*Phi(:).');% even azimuth
    % Field intensity (Richards–Wolf combination):
    % See eg. Sheppard/Török: Ix + Iy from transverse + longitudinal
    % A common compact form is:
    Ex2 = abs(I0 + I2).^2;
    Ey2 = abs(I0 - I2).^2;
    Ez2 = abs(2*I1).^2;
    I   = real(Ex2 + Ey2 + Ez2);         % 1×P

    PSF(:,:,iz) = reshape(I, size(rho));
end

% downsample to camera pixels
if params.sf>1
    PSF = bin2(PSF, params.sf);
end

if wantDXP
    % Numerical derivative w.r.t xp using symmetric difference:
    h = max(dx_cam/4, 1e-10);
    outL = vectorialPSF(xp-h, zvec, nx, params);
    outR = vectorialPSF(xp+h, zvec, nx, params);
    dPSFdxp = (outR.PSF - outL.PSF) / (2*h);
end

out.PSF  = PSF;
if wantDXP, out.dxp = dPSFdxp; end
out.meta = params;
end

% -------------------- helpers --------------------
function [x,w]=gaussLegendre(a,b,n)
% nodes/weights on [a,b]
% uses poly recursion (stable for n<=1000)
beta = 0.5 ./ sqrt(1 - (2*(1:n-1)).^-2);
T = diag(beta,1)+diag(beta,-1);
[V,D] = eig(T);
x = (b-a)/2*diag(D)'+(a+b)/2;
w = (b-a)/2*(2*(V(1,:).^2));
[x,ix]=sort(x); w=w(ix);
end

function v=getfielddef(s, f, dv)
if isfield(s,f), v = s.(f); else, v = dv; end
end

function Y = bin2(X, sf)
% average non-overlapping sf×sf blocks in x,y
[m,n,k] = size(X);
M = floor(m/sf)*sf; N=floor(n/sf)*sf;
X = X(1:M,1:N,:);
Y = squeeze(mean(reshape(X,M/sf,sf,N/sf,sf,k),[2 4]));
end
