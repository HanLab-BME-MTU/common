function [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi, delta_x, delta_y, i_end, j_end)
% LSFIRSTSECONDDIFFERENCES calculates gradient wit hsecond order accuracy
%    
%
%
% SYNOPSIS   [delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi, delta_x, delta_y, i_end, j_end)
%
%
% INPUT      phi        : phi=f(x,y) function values on a grid
%            delta_x    : x-direction grid spacing
%            delta_y    : y-direction grid spacing
%            i_end      : number of x grid points
%            j_end      : number of y grid points 
%                          
% 
% OUTPUT     delta_plus     :  right side absolute value of the gradient
%            delta_minus    :  left side absolute value of the gradient
%            grad_x         :  x-comp. of the gradient
%            grad_y         :  y-comp. of the gradient
%                           
% DEPENDENCES     lsGradient2o uses {                                
%                                       }
%
%                 lsGradient2o is used by { 
%                                           }
%
% Matthias Machacek 06/22/04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%% Difference operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dy_minus          =  zeros(i_end, j_end);
Dy_plus           =  zeros(i_end, j_end);
Dy_minus_y_minus  =  zeros(i_end, j_end);
Dy_plus_y_plus    =  zeros(i_end, j_end);
Dy_plus_y_minus   =  zeros(i_end, j_end);

Dx_minus          =  zeros(i_end, j_end);
Dx_plus           =  zeros(i_end, j_end);
Dx_minus_x_minus  =  zeros(i_end, j_end);
Dx_plus_x_plus    =  zeros(i_end, j_end);
Dx_plus_x_minus   =  zeros(i_end, j_end);

delta_x_sq = delta_x^2;
delta_y_sq = delta_y^2;

for i=1:i_end
   for j=1:j_end
      if i < 3
         % assume that the value at i=1 is the same as at i=2
         Dy_minus(i,j)          =  (phi(i+1,j)-  phi(i  ,j)           )  / delta_y;
         Dy_plus(i,j)           =  (phi(i+1,j)-  phi(i  ,j)           )  / delta_y;
         Dy_minus_y_minus(i,j)  =  (phi(i+2,j)-2*phi(i+1,j)+phi(i,j)  )  / delta_y_sq;
         Dy_plus_y_plus(i,j)    =  (phi(i+2,j)-2*phi(i+1,j)+phi(i,j)  )  / delta_y_sq;
         Dy_plus_y_minus(i,j)   =  (phi(i+2,j)-2*phi(i+1,j)+phi(i,j)  )  / delta_y_sq;         
      elseif i > i_end-3
         Dy_minus(i,j)          =  (phi(i  ,j)-  phi(i-1,j)             )  / delta_y;
         Dy_plus(i,j)           =  (phi(i  ,j)-  phi(i-1,j)             )  / delta_y;
         Dy_minus_y_minus(i,j)  =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_y_sq;
         Dy_plus_y_plus(i,j)    =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_y_sq;
         Dy_plus_y_minus(i,j)   =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_y_sq;       
      else
         Dy_minus(i,j)          =  (phi(i  ,j)-  phi(i-1,j)             )  / delta_y;
         Dy_plus(i,j)           =  (phi(i+1,j)-  phi(i  ,j)             )  / delta_y;
         Dy_minus_y_minus(i,j)  =  (phi(i  ,j)-2*phi(i-1,j)+phi(i-2,j)  )  / delta_y_sq;
         Dy_plus_y_plus(i,j)    =  (phi(i+2,j)-2*phi(i+1,j)+phi(i  ,j)  )  / delta_y_sq;
         Dy_plus_y_minus(i,j)   =  (phi(i+1,j)-2*phi(i  ,j)+phi(i-1,j)  )  / delta_y_sq;      
      end
      
      if j < 3
         Dx_minus(i,j)          =  (phi(i,j+1)-  phi(i,j  )             )  / delta_x;
         Dx_plus(i,j)           =  (phi(i,j+1)-  phi(i,j  )             )  / delta_x;
         Dx_minus_x_minus(i,j)  =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_x_sq;
         Dx_plus_x_plus(i,j)    =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_x_sq;
         Dx_plus_x_minus(i,j)   =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_x_sq;        
      elseif j > j_end-3
         Dx_minus(i,j)          =  (phi(i,j  )-  phi(i,j-1)             )  / delta_x;
         Dx_plus(i,j)           =  (phi(i,j  )-  phi(i,j-1)             )  / delta_x;
         Dx_minus_x_minus(i,j)  =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_x_sq;
         Dx_plus_x_plus(i,j)    =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_x_sq;
         Dx_plus_x_minus(i,j)   =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_x_sq;    
      else
         Dx_minus(i,j)          =  (phi(i,j  )-  phi(i,j-1)             )  / delta_x;
         Dx_plus(i,j)           =  (phi(i,j+1)-  phi(i,j  )             )  / delta_x;
         Dx_minus_x_minus(i,j)  =  (phi(i,j  )-2*phi(i,j-1)+phi(i,j-2)  )  / delta_x_sq;
         Dx_plus_x_plus(i,j)    =  (phi(i,j+2)-2*phi(i,j+1)+phi(i,j  )  )  / delta_x_sq;
         Dx_plus_x_minus(i,j)   =  (phi(i,j+1)-2*phi(i,j  )+phi(i,j-1)  )  / delta_x_sq;
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

grad_x =        zeros(i_end, j_end);
grad_y =        zeros(i_end, j_end);
delta_plus  =   zeros(i_end, j_end);
delta_minus =   zeros(i_end, j_end);

delta_x_half = delta_x/2;
delta_y_half = delta_y/2;

for i=1:i_end
    for j=1:j_end
%         if  (Dx_minus_x_minus(i,j) *  Dx_plus_x_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dx_minus_x_minus(i,j)) <= abs(Dx_plus_x_minus(i,j) < 0)
%             m = Dx_minus_x_minus(i,j);
%         else
%             m = Dx_plus_x_minus(i,j);
%         end
%         A = Dx_minus(i,j) + delta_x_half * m;
%         
%         if  (Dx_plus_x_plus(i,j) * Dx_plus_x_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dx_plus_x_plus(i,j)) <= abs(Dx_plus_x_minus(i,j))
%             m = Dx_plus_x_plus(i,j);
%         else
%             m = Dx_plus_x_minus(i,j);
%         end
%         B = Dx_plus(i,j)  - delta_x_half * m;
%         
%         if  (Dy_minus_y_minus(i,j) *  Dy_plus_y_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dy_minus_y_minus(i,j)) <= abs(Dy_plus_y_minus(i,j))
%             m = Dy_minus_y_minus(i,j);
%         else
%             m = Dy_plus_y_minus(i,j);
%         end      
%         C = Dy_minus(i,j) + delta_y_half * m;
%         
%         if  (Dy_plus_y_plus(i,j) * Dy_plus_y_minus(i,j)) < 0
%             m = 0;
%         elseif abs(Dy_plus_y_plus(i,j)) <= abs(Dy_plus_y_minus(i,j))
%             m = Dy_plus_y_plus(i,j);
%         else
%             m = Dy_plus_y_minus(i,j);
%         end         
%         D = Dy_plus(i,j)  - delta_y_half * m;
       
       
       
      A = Dx_minus(i,j) + delta_x_half * switch_m(Dx_minus_x_minus(i,j), Dx_plus_x_minus(i,j));
      B = Dx_plus(i,j)  - delta_x_half * switch_m(Dx_plus_x_plus(i,j),   Dx_plus_x_minus(i,j));
      C = Dy_minus(i,j) + delta_y_half * switch_m(Dy_minus_y_minus(i,j), Dy_plus_y_minus(i,j));
      D = Dy_plus(i,j)  - delta_y_half * switch_m(Dy_plus_y_plus(i,j),   Dy_plus_y_minus(i,j));
      
      grad_x(i,j) = max(A,0) + min(B,0);
      grad_y(i,j) = max(C,0) + min(D,0);
      
      delta_plus(i,j)  =  sqrt(max(A,0)^2 + min(B,0)^2+...
                               max(C,0)^2 + min(D,0)^2);
      
      delta_minus(i,j) =  sqrt(max(B,0)^2 + min(A,0)^2+...
                               max(D,0)^2 + min(C,0)^2);      
   end
end



function m = switch_m(D1, D2)

if (D1 * D2) < 0
   m = 0; 
elseif abs(D1) <= abs(D2)
   m = D1;
else
   m = D2;
end