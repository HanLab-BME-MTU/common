p.lambda = 520e-9;
p.ni = 1.515; p.ng = 1.525; p.ns = 1.333;
p.NA = 1.4; p.M = 100; p.pixel = 6.5e-6;
p.t_ig = 0; p.t_gs = 170e-6;      % example coverslip thickness
p.Ntheta = 400; p.sf = 3;

nx   = 129;
zvec = linspace(-500e-9, 500e-9, 41);
xp   = 0;

out = vectorialPSFnew(xp, zvec, nx, p);
figure,
imagesc(out.PSF(:,:,ceil(end/2))); axis image off; colormap hot
title('PSF at focus (z≈0)')