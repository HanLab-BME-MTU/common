function F = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x, grad_y, domain)
% LSGETVELOCITYFCT calculates the speed F at the grid points
%    
%
%
% SYNOPSIS   F = lsGetVelocityFct(phi, phi_target, i_end, j_end, grad_x, grad_y, domain)
%
%
% INPUT      phi        :
%            phi_target :
%            i_end      :
%            j_end      :
%            grad_x     :
%            grad_y     :
%            domain     :
%                          
% 
% OUTPUT     F         :
%              
%                           
% DEPENDENCES    lsGetVelocityFct uses {                                
%                                       }
%
%                lsGetVelocityFct is used by { 
%                                           }
%
% Matthias Machacek 06/24/04

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
         grad_l = sqrt(grad_x(i,j)^2 + grad_y(i,j)^2);
         if grad_l == 0
             grad_l = 0.000001;
         end
         
         % x-translation
         %F(i,j) =  grad_x(i,j) / grad_l;
         
         % y-translation
         % F(i,j) =  grad_y(i,j) / grad_l;
         
         % rotation
         F(i,j) =  -j/2*grad_x(i,j) / grad_l + i/2*grad_y(i,j) / grad_l;
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
elseif 1
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

    % level set difference driven flow 
    kappa = lsCurvature(phi, domain.x_spacing, domain.y_spacing, i_end, j_end);
    for i=1:i_end
        for j=1:j_end
            %F(i,j) = - sign(phi_target(i,j))*2;
            
            %sign_vel = sign(phi(i,j) - phi_target(i,j));
            %F(i,j) = sign_vel *( phi(i,j) - phi_target(i,j))^2;
            %F(i,j) = sign( phi(i,j) - phi_target(i,j)) * log(phi(i,j) - phi_target(i,j));
            % this works
            %F(i,j) = atan(phi(i,j) - phi_target(i,j));
            d_level = phi(i,j) - phi_target(i,j);
            if d_level >= 0
                % protrusion
                F(i,j) = 0.05 * d_level + atan(d_level);
            else
                % retraction
                F(i,j) = d_level * (0.5 + kappa(i,j));
            end
        end
    end
else
    % Curvature driven flow
    F = lsCurvature(phi, domain.x_spacing, domain.y_spacing, i_end, j_end);
    %F = 1 - 0.35 .* F;
    F = - F;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%