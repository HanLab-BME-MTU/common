% Synthetic 2-D Gaussian
nx = 21; ny = 21;
[xg, yg] = meshgrid( (-floor(nx/2):floor(nx/2)), (-floor(ny/2):floor(ny/2)) );
true_xp = 0.7; true_yp = -1.2;
true_A  = 540;
true_s  = 1.8;
true_c  = 20;

G = true_A*exp(-((xg-true_xp).^2 + (yg-true_yp).^2)/(2*true_s^2)) + true_c;
G = G + randn(size(G))*2.0;  % a touch of noise

% Initial guess & mode (estimate xp,yp,A,s,c)
p0 = [0 0 max(G(:)) 2 min(G(:))];
mode = 'xyasc';                      % order: [xp yp A s c]
opts = [200 1e-6 1e-6];              % [maxIter eAbs eRel]

% Call your wrapper
[prm, prmStd, C, res, J] = fitGaussian2D(G, p0, mode, opts);

disp('Estimated [xp yp A s c] vs true:');
disp(prm);
disp([true_xp true_yp true_A true_s true_c]);

% Basic checks
assert(numel(prm)==5);
assert(numel(prmStd)==sum(ismember(mode, 'xyasc')));
assert(isequal(size(J), [numel(G(~isnan(G))) sum(ismember(mode,'xyasc'))]));
fprintf('RSS: %.3f  mean(res): %.3f  std(res): %.3f\n', res.RSS, res.mean, res.std);