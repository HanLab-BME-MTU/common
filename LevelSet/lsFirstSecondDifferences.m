function [delta_plus, delta_minus, grad_x, grad_y] = lsFirstSecondDifferences(phi, delta_x, delta_y, i_end, j_end)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%% Difference operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:i_end
   for j=1:j_end
      if i < 3
         % assume that the value at i=1 is the same as at i=2
         Dx_minus(i,j)          =  (phi(i+1,j)-  phi(i  ,j)           )  / delta_x;
         Dx_plus(i,j)           =  (phi(i+1,j)-  phi(i  ,j)           )  / delta_x;
         Dx_minus_x_minus(i,j)  =  (phi(i+2,j)-2*phi(i+1,j)+phi(i,j)  )  / delta_x^2;
         Dx_plus_x_plus(i,j)    =  (phi(i+2,j)-2*phi(i+1,j)+phi(i,j)  )  / delta_x^2;
         Dx_plus_x_minus(i,j)   =  (phi(i+2,j)-2*phi(i+1,j)+phi(i,j)  )  / delta_x^2;         
      elseif i > i_end-3
         Dx_minus(i,j)          =  (phi(i  ,j)-  phi(i-1,j)             )  / delta_x;
         Dx_plus(i,j)           =  (phi(i  ,j)-  phi(i-1,j)             )  / delta_x;
         Dx_minus_x_minus(i,j)  =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_x^2;
         Dx_plus_x_plus(i,j)    =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_x^2;
         Dx_plus_x_minus(i,j)   =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_x^2;       
      else
         Dx_minus(i,j)          =  (phi(i  ,j)-  phi(i-1,j)             )  / delta_x;
         Dx_plus(i,j)           =  (phi(i+1,j)-  phi(i  ,j)             )  / delta_x;
         Dx_minus_x_minus(i,j)  =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_x^2;
         Dx_plus_x_plus(i,j)    =  (phi(i+2,j)-2*phi(i+1,j)+phi(i  ,j)  )  / delta_x^2;
         Dx_plus_x_minus(i,j)   =  (phi(i+1,j)-2*phi(i  ,j)+phi(i-1,j)  )  / delta_x^2;      
      end
      
      if j < 3
         Dy_minus(i,j)          =  (phi(i,j+1)-  phi(i,j  )             )  / delta_y;
         Dy_plus(i,j)           =  (phi(i,j+1)-  phi(i,j  )             )  / delta_y;
         Dy_minus_y_minus(i,j)  =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_y^2;
         Dy_plus_y_plus(i,j)    =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_y^2;
         Dy_plus_y_minus(i,j)   =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_y^2;        
      elseif j > j_end-3
         Dy_minus(i,j)          =  (phi(i,j  )-  phi(i,j-1)             )  / delta_y;
         Dy_plus(i,j)           =  (phi(i,j  )-  phi(i,j-1)             )  / delta_y;
         Dy_minus_y_minus(i,j)  =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_y^2;
         Dy_plus_y_plus(i,j)    =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_y^2;
         Dy_plus_y_minus(i,j)   =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_y^2;    
      else
         Dy_minus(i,j)          =  (phi(i,j  )-  phi(i,j-1)             )  / delta_y;
         Dy_plus(i,j)           =  (phi(i,j+1)-  phi(i,j  )             )  / delta_y;
         Dy_minus_y_minus(i,j)  =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_y^2;
         Dy_plus_y_plus(i,j)    =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_y^2;
         Dy_plus_y_minus(i,j)   =  (phi(i,j+1)-2*phi(i,j  )+phi(i,j-1)  )  / delta_y^2;
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


for i=1:i_end
   for j=1:j_end
      A(i,j) = Dx_minus(i,j) + delta_x/2 * switch_m(Dx_minus_x_minus(i,j), Dx_plus_x_minus(i,j));
      B(i,j) = Dx_plus(i,j)  - delta_x/2 * switch_m(Dx_plus_x_plus(i,j),   Dx_plus_x_minus(i,j));
      C(i,j) = Dy_minus(i,j) + delta_y/2 * switch_m(Dy_minus_y_minus(i,j), Dy_plus_y_minus(i,j));
      D(i,j) = Dy_plus(i,j)  - delta_y/2 * switch_m(Dy_plus_y_plus(i,j),   Dy_plus_y_minus(i,j));
      
      grad_x(i,j) = max(A(i,j),0) + min(B(i,j),0);
      grad_y(i,j) = max(C(i,j),0) + min(D(i,j),0);
      
      delta_plus(i,j)  =  sqrt(max(A(i,j),0)^2 + min(B(i,j),0)^2+...
                               max(C(i,j),0)^2 + min(D(i,j),0)^2);
      
      delta_minus(i,j) = sqrt(max(B(i,j),0)^2 + min(A(i,j),0)^2+...
                              max(D(i,j),0)^2 + min(C(i,j),0)^2);      
   end
end



function m = switch_m(D1, D2)

if (D1 * D2) < 0
   m=0; 
elseif abs(D1) <= abs(D2)
   m = D1;
else
   m = D2;
end