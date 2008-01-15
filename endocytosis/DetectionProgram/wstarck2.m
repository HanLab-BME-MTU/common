function y=wstarck2(x,k,irecon)

% Reconstructs the image without noise by summing all significant
% coefficients in each detail

[ly,lx]=size(x);		% Get the original image size


%------------------------------
%   SETUP ARRAYS 
%------------------------------

y=zeros(ly,lx);
m=zeros(ly,lx,k);


tx=x;				% Copy the original.
yl(:,:,1)=x;
w=zeros(ly,lx);

%------------------------------
%     START THE ALGORITHM 
%------------------------------
s=zeros(ly,lx);


for ind=1:k			% For every scale...
    
    
    
    %%% Now let's call wtlo2 (low pass transform) first over the rows %%%
    %%% of the present image			   %%%
    
    
    tw(1:ly,1:lx)=wtlo2(tx(1:ly,1:lx),ind);
    
    
    %%% And then over the columns  %%%
    
    
    tw(1:ly,1:lx)=wtlo2(tw(1:ly,1:lx)',ind)';
    
    
    yl(:,:,ind+1)=tw;	% This is the ind_th approximation of x
    
    
    yh(:,:)=yl(:,:,ind)-yl(:,:,ind+1); % This is the ind_th detail or residue
    %yh=tx-tw; % This is the ind_th detail or residue (or high pass component)
    tx=tw;   % This is the ind_th approximation of x (or low pass component)   
    st=std(yh(:)); % Estimate the standard deviation of the noise in the detail at each scale.
    for i=1:ly
        for j=1:lx
            %if abs(yh(i,j)) >=(4-floor((ind+1)/2))*st
            if abs(yh(i,j)) >=3*st
                s(i,j)=1; %If s=0 for all scales then the pixel has mostly noise
            end
        end
    end
    w=w+s.*yh;
end				% end of all scales.
if irecon==1
    w=w+tx; %filtered image with last approximation added
    %elseif irecon==2
    %    w=s.*x; % multiresolution transform
end
%ttmn=min(min(w))
%ttmx=max(max(w))

tx=x-w; %error image

sgm2=std(tx(:));


y=w;  

%------------------------------
%    END OF THE ALGORITHM 
%------------------------------
