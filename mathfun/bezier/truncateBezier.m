function [cP,t] = truncateBezier(cP,tStart,tEnd,t)
% function [cP,t] = truncateBezier(cP,tStart,tEnd,t)
% truncateBezier computes the control points of a segment of the input 
% Bezier curve. The input curve can be a linear, quadratic or cubic 
% Bezier curve. Optionally, points of the old curve can be mapped onto the
% new curve.
%
% Required Inputs:
% cP                A N x D array representing a set of 3-dimensional 
%                   control points where N is 2,3 or 4. D is the dimension
%                   of the control points.
% tStart            Defines the start point of the segment            
% tEnd              Defines the end point of the segment
% 
% Optional Inputs:
% t                 M x 1 array containing parameter values of the old
%                   curve that will be mapped to the new one
%
% Outputs:
% cP                The control points of the new curve
% t                 The transformed parametrization 
%
% Pascal Berard, June 2011

% Transform the parametrization
if nargin == 4
    t = (t-tStart)/(tEnd-tStart);
end

% Compute the new control points
switch (size(cP,1))
    case 2
        % Linear Bezier
        % B[t_] = (1 - t)*P0 + t*P1;
        % B[t_] = a*t + b;
        a = -cP(1,:) + cP(2,:);
        b = cP(1,:);
        
        % P0T -> b + a tStart
        cP(1,:) = b + a*tStart;
        
        % P1T -> b + a tEnd
        cP(2,:) = b + a*tEnd;
        
    case 3
        % Quadratic Bezier
        % B[t_] = (1 - t)^2*P0 + 2*(1 - t)*t*P1 + t^2*P2;
        % B[t_] = a*t^2 + b*t + c;
        a = cP(1,:) - 2*cP(2,:) + cP(3,:);
        b = -2*cP(1,:) + 2*cP(2,:);
        c = cP(1,:);
        
        % P0T -> c + b tStart + a tStart^2
        cP(1,:) = c + b*tStart + a*tStart^2;
        
        % P1T -> c + (b tEnd)/2 + (b tStart)/2 + a tEnd tStart
        cP(2,:) = c + (b*tEnd)/2 + (b*tStart)/2 + a*tEnd*tStart;
        
        % P2T -> c + b tEnd + a tEnd^2
        cP(3,:) = c + b*tEnd + a*tEnd^2;
        
    case 4
        % Cubic Bezier
        % B[t_] = (1 - t)^3*P0 + 3*(1 - t)^2*t*P1 + 3*(1 - t)*t^2*P2 + t^3*P3;
        % B[t_] = a*t^3 + b*t^2 + c*t + d;
        a = -cP(1,:) + 3*cP(2,:) - 3*cP(3,:) + cP(4,:);
        b = 3*cP(1,:) - 6*cP(2,:) + 3*cP(3,:);
        c = -3*cP(1,:) + 3*cP(2,:);
        d = cP(1,:);
        
        % P0T -> d + c tStart + b tStart^2 + a tStart^3
        cP(1,:) = d + c*tStart + b*tStart^2 + a*tStart^3;
        
        % P1T -> d + (c tEnd)/3 + (2 c tStart)/3 + (2 b tEnd tStart)/3 + (b tStart^2)/3 + a tEnd tStart^2
        cP(2,:) = d + (c*tEnd)/3 + (2*c*tStart)/3 + (2*b*tEnd*tStart)/3 + (b*tStart^2)/3 + a*tEnd*tStart^2;
        
        % P2T -> d + (2 c tEnd)/3 + (b tEnd^2)/3 + (c tStart)/3 + (2 b tEnd tStart)/3 + a tEnd^2 tStart
        cP(3,:) = d + (2*c*tEnd)/3 + (b*tEnd^2)/3 + (c*tStart)/3 + (2*b*tEnd*tStart)/3 + a*tEnd^2*tStart;
        
        % P3T -> d + c tEnd + b tEnd^2 + a tEnd^3
        cP(4,:) = d + c*tEnd + b*tEnd^2 + a*tEnd^3; 
    otherwise
        assert(false,'Only linear, quadratic, and cubic Bézier 3D curves supported!');
end

end

