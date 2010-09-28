function wm=wstarck222(x,k,irecon)

% Calculates the multiscale product from the reconstructed image

[ly,lx]=size(x);		% Get the original image size


%------------------------------
%   SETUP ARRAYS
%------------------------------
yl = zeros(ly,lx,k);

tx=x;		% Copy the original.
yl(:,:,1)=x;
wm = ones(ly,lx);

for ind = 1:k			% For every scale...
    
    %%% Now let's call wtlo2 (low pass transform) first over the rows %%%
    %%% of the present image			   %%%
    yl(:,:,ind+1) = wtlo2(wtlo2(tx, ind)',ind)';	% This is the ind_th approximation of x
    
    %wi = awt(tx, ind);
    %yl(:,:,ind+1) = wi(:,:,end);
    
    
    yh = yl(:,:,ind) - yl(:,:,ind+1); % This is the ind_th detail or residue
    wm = abs(wm.*yh);
    %yh=tx-yl(:,:,ind+1); % This is the ind_th detail or residue (or high pass component)
    tx = yl(:,:,ind+1);   % This is the ind_th approximation of x (or low pass component)
    %st=std(yh(:)); % Estimate the standard deviation of the noise in the detail at each scale.
end				% end of all scales.

