function w = wstarck2(x,k,irecon)
% WSTARCK2 reconstructs the image without noise by summing all significant
% coefficients in each detail.  
%  Y = WSTARCK2(X,K,IRECON), where X is the image to be processed.
%    K is the level of decomposition for wavelet transformation.
%    IRECON is the flag to indicate if the ind_th approximation of X (or
%    low pass component) should be added.  When IRECON = 0, WSTARCK2 just
%    calcualted reconstructed image up to order k without kth detail.
%
%    Reference:  "Olivo-Marin J.C. 2002. Extraction of spots in biological
%    images using multiscale products. Pattern Recognit. 35: 1989�1996."
% 
% Last updated:  Shann-Ching Sam Chen, April 8, 2008
%                add description to the function and rewrote the inner loop 
%                in a vectorized form to speed up computation

[ly,lx]=size(x);		% Get the original image size

%------------------------------
%   SETUP ARRAYS
%------------------------------
yl = zeros(ly, lx, k + 1);

tx = x;
yl(:,:,1) = x;
w = zeros(ly,lx);
s = zeros(ly,lx);

for ind = 1:k
    yl(:,:,ind+1) = wtlo2(wtlo2(tx, ind)', ind)';
    %wi = awt(tx, ind);
    %yl(:,:,ind+1) = wi(:,:,end);
    
    yh = yl(:,:,ind) - yl(:,:,ind+1); % This is the ind_th detail or residue
    tx = yl(:,:,ind+1);   % This is the ind_th approximation of x (or low pass component)  
    s(abs(yh)>=3*std(yh(:))) = 1;
    
    w = w + s.*yh;
end

if irecon==1
    w = w + tx;
end