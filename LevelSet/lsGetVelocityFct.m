function F = lsGetSpeedFct(phi, phi_target, i_end, j_end, grad_x, grad_y, domain)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculate the speed function F  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
   % Uniform velocity field
   for i=1:i_end
      for j=1:j_end
         F(i,j) = 1;
      end
   end
elseif 0
   % Flow in x-direction, y- direction, circular velocity field
   for i=1:i_end
      for j=1:j_end
         grad_phi_l = sqrt(grad_x^2 + grad_y^2);
         
         F(i,j) =  -j/10*grad_phi_x / grad_phi_l + i/10*grad_phi_y / grad_phi_l;
      end
   end
elseif 0
   % Gradient driven flow
   for i=1:i_end
      for j=1:j_end
         grad_phi_x = max(A(i,j),0) + min(B(i,j),0);
         grad_phi_y = max(C(i,j),0) + min(D(i,j),0);
         grad_phi_l = sqrt(grad_phi_x^2 + grad_phi_y^2);
         
         F(i,j) =  -Dx_plus_target(i,j) * grad_x - Dy_plus_target(i,j) * grad_y;
      end
   end  
else
   % Distance driven flow
   if 0
   % Calculate the velocity at the front, thus at the zero level set
   % find the zero level set:
   phi_zero = lsGetZeroLevel(phi, domain);
   % get the velocity at the zero level set
   for i = 1:size(phi_zero,2)
       F(phi_zero(1,i), phi_zero(2,i)) = 0;
   end
   
   % get the extended velocity
   end
   
   for i=1:i_end
      for j=1:j_end
         F(i,j) = - sign(phi_target(i,j))*2;
      end
   end    
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%