function [phi_next, F, delta_t_optimal] = lsSolveConvection(phi, delta_t, delta_x, delta_y, i_end, j_end, phi_target, domain)
% LSSOLVECONVECTION solves the equation d phi/ dt + F|grad phi| = 0 
%    
%           Second order spatial
%           and first order temporal scheme 
%           according to J.A. Sethian: Level Set Methods
%           and Fast Marching Methods p.: 66
%
%
% SYNOPSIS   [phi_next, velocity_fct, delta_t_optimal] = lsSolveConvection(phi, delta_t, delta_x, delta_y, i_end, j_end, phi_target, domain)
%
%
% INPUT      phi        :
%            delta_t    :
%            delta_x    :
%            delta_y    :
%            i_end      :
%            j_end      :
%            phi_target :
%            domain     :
%                          
% 
% OUTPUT     phi_next       :
%            velocity_fct   : used velocity fct
%            delta_t_optimal: optimal time step 
%              
%                           
% DEPENDENCES    lsSolveConvection uses {                                
%                                       }
%
%                lsSolveConvection is used by { 
%                                           }
%
% Matthias Machacek 06/22/04


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%% Difference operators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[delta_plus, delta_minus, grad_x, grad_y] = lsGradient2o(phi, delta_x, delta_y, i_end, j_end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Difference operators for the target function %%%%%%
for i=1:i_end
   for j=1:j_end
      if i > i_end-1  
         Dx_plus_target(i,j)           =  (phi_target(i  ,j)-phi_target(i-1,j)) / delta_x;  
      else
         Dx_plus_target(i,j)           =  (phi_target(i+1,j)-phi_target(i  ,j)) / delta_x;
      end
      
      if j > j_end-1
         Dy_plus_target(i,j)           =  (phi_target(i,j  )-phi_target(i,j-1)) / delta_y;
      else
         Dy_plus_target(i,j)           =  (phi_target(i,j+1)-phi_target(i,j  )) / delta_y;
      end
   end
end

% get the gradient at the zero level set position of the target function
for i=1:i_end
   for j=1:j_end
      grad_phi_target(i,j) = sqrt(Dx_plus_target(i,j)^2+Dy_plus_target(i,j)^2);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculate the speed function F  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x, grad_y, domain);


% get the maximal time step based on the CLF number 
delta_t_optimal = 0.9 * delta_x / max(max(F));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the level-sets   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% methods Euler, Adams, Heun

method = 'Heun';  

if strcmp(method,'Euler')
   % First order acccurate: Euler method
   for i=1:i_end
      for j=1:j_end
         phi_next(i,j) = phi(i,j) - delta_t*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j));
      end
   end
elseif strcmp(method,'Adams')
   % Adams_Bashforth method: second order
    for i=1:i_end
      for j=1:j_end
         phi_next(i,j) = phi(i,j) - delta_t*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j));
      end
   end  
elseif strcmp(method,'Heun')
   % Heun method: second order. A predictor-corrector schema
   % Predictor step
   for i=1:i_end
      for j=1:j_end
         phi_p(i,j) = phi(i,j) - delta_t*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j));
      end
   end
   
   [delta_plus_p, delta_minus_p, grad_x_p, grad_y_p] = lsFirstSecondDifferences(phi_p,...
                                                                     delta_x, delta_y, i_end, j_end);
   
   F_p = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x_p, grad_y_p, domain);
   
   % Correctror step
   for i=1:i_end
      for j=1:j_end
         phi_next(i,j) = phi(i,j) - delta_t/2*(max(F(i,j), 0)*delta_plus(i,j) + min(F(i,j), 0) * delta_minus(i,j)+...
                                               max(F_p(i,j), 0)*delta_plus_p(i,j) + min(F_p(i,j), 0) * delta_minus_p(i,j));
      end
   end   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%