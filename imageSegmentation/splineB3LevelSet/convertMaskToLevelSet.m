function phi = convertMaskToLevelSet(mask)

e = edge(mask);
phi = double(bwdist(double(e)));
phi(mask == false) = -phi(mask == false);

% % TODO: implement RK2 or RK3 if needed
% % TODO: use higher order numerical scheme for spatial gradient approximate
% % TODO: implement the volume constraint (sussman 1998)
% % TODO: speed up computation (implement in C++)
% 
% phi0 = zeros(size(mask));
% 
% indIn = find(mask == true);
% indOut = find(mask == false);
% 
% phi0(indIn) = 1;
% phi0(indOut) = -1;
% 
% G = zeros(size(mask));
% 
% iter = 0;
% 
% phi = phi0;
% 
% zIn = zeros(size(indIn));
% zOut = zeros(size(indOut));
% 
% dt = .5;
% 
% maxIter = max(size(mask));
% 
% while iter < maxIter
%     
%     phiPadded = padarray(phi, [1 1], 'replicate');
%     
%     A = phi - phiPadded(2:end-1,1:end-2);
%     B = phiPadded(2:end-1,3:end) - phi;
%     C = phi - phiPadded(1:end-2,2:end-1);
%     D = phiPadded(3:end,2:end-1) - phi;
%     
%     G(indIn) = sqrt(max(max(A(indIn),zIn).^2, min(B(indIn),zIn).^2) + ...
%         max(max(C(indIn),zIn).^2, min(D(indIn),zIn).^2)) - 1;
%     
%     G(indOut) = sqrt(max(min(A(indOut),zOut).^2, max(B(indOut),zOut).^2) + ...
%         max(min(C(indOut),zOut).^2, max(D(indOut),zOut).^2)) - 1;
%     
%     phi = phi - dt * phi0 .* G;
%     
%     iter = iter + 1;
% end
