function y=wtlo2(x,k)

[liy,lix]=size(x);

%------------------------------
%     START THE ALGORITHM 
%------------------------------

for it=1:liy		% For every row of the input matrix...
			% (this makes one wavelet transform
			% for each of the rows of the input matrix)
    edge=2^(k+1);
    lmax=lix+2*edge;
    % Copy the vector to transform.
    t =zeros(1,lmax);
     
   
	for jx=1:lmax   % apply continuous boundary conditions
        if jx < edge+1
            %t(jx)=x(it,lix-edge+jx);
            t(jx)=x(it,1);
        elseif jx > lix+edge
            %t(jx)=x(it,jx-lix-edge);
            t(jx)=x(it,lix);
        else 
            t(jx)=x(it,jx-edge);
        end
    end
                 
    %t=x(it,:);			% Copy the vector to transform.
    
        
        for ix=1+edge:lix+edge
            % Do lowpass filtering ...
            % first line below gives the better B3-spline a trous filter
            yl(ix-edge)=0.375*t(ix)+0.25*t(ix+2^(k-1))+0.25*t(ix-2^(k-1))+0.0625*t(ix+2^k)+0.0625*t(ix-2^k);
            % next line gives simplest (triangle) a trous filter
            %yl(ix-edge)=0.5*t(ix)+0.25*t(ix+2^(k-1))+0.25*t(ix-2^(k-1));
                      
        end
        
    	
	%end

	y(it,:)=yl;		       	% Wavelet vector (1 row vector)
    
end				% End of the "rows" loop.

%------------------------------
%    END OF THE ALGORITHM 
%------------------------------