function interField=vecFldInterpAnisoA(M, tgtPoints, flowVecList, sigmaU, sigmaV, polygon)
% vecFldInterpAnisoA uses a user-defined directional anisotropic filter to interpolate
% the vector field on a user-specified point set based on flow vectors on a given point set.
% Directions of the anisotropic filters are provided externally. This function will be
% useful if the global flow field is already known.
%
% If no external flow field information is directly available, another related
% function "vecFldInterpAnisoB" can be used. This function will
% itereatively calculate the local flow field based on given flow vectors.
%
% Interpolation is based on an anisotropic Gaussain probability density function of the following form:
%  
% p(u, v) =  1 / (2 * pi * sigmaU * sigmaV) * exp(-0.5 * (u / sigmaU)^2 -
%            0.5 * (v / sigmaV)^2);
%
% SYNOPSIS   interpFiled = vecFldInterpAnisoA(M, tgtPoints, flowVecList, sigmaU, sigmaV, R, polygon)
%
% INPUT              M  :     known vector field, stored in a (nx4)-matrix of the form [y0 x0 y x]n
%                             (where (y0,x0) is the base and (y,x) is the tip of
%                             the vector).
%
%            tgtPoints  :     points at which vector field needs to be calculated.
%                             It is stored in a (mx2)-matrix of the form
%                             [yg xg]m. In particular, this point set can be
%                             the same as the given point set in M.
%
%          flowVecList  :     list of flow vectors for the corresponding
%                             points in tgtPoints. These vectors are not necessarily
%                             normalized.
%             
%       sigmaU, sigmaV  :     Gauss distribution parameters for calculating
%                             the probability function. Currently, it is
%                             assumed that they are contant over the entire
%                             field. sigmaU is defined along the flow
%                             direction. sigmaV is defined perpendicular to
%                             the flow direction.
%                             NOTICE: these two parameters will also be used to determine a search region, based 
%                             essentially on the three-sigma rule. 
% 
%              polygon  :     (optional - pass polygon=[] to disable). The interpolated vector
%                             can be cropped to remove vectors outside a given region of interest.
%                             To create the polygon use the functions ROIPOLY or
%                             GETLINE. These functions return the polygon vertices
%                             stored in two vectors y and x. Set polygon=[y x] to
%                             use with vectorFieldInterp.
%
% OUTPUT    interField  :     interpolated vector field for the points given in tgtPoints.

% DEPENDENCES                 This function calls vecFldInterpAnisoA
%
% REMARKS                     This implementation is based on the function and
%                             interface of "vectorFieldInterp" created by Aaron Ponti,
%                             11/18/2002. 
%
% AUTHOR                      Ge Yang, Ph.D. 
%                             geyang@scripps.edu
%                             Laboratory for Computational Cell Biology
%                             The Scripps Research Institute
%
% DATE OF CREATION            March 10, 2004
% 
% CHANGES
%                             ver 0.3: The code is reorganized for better
%                             readability and initial release.
%                               
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% Vector base positions
srcPoints = M(:, 1:2); % These will be the centers of the probability ellipsoids.
V=[M(:, 3) - M(:, 1) M(:, 4) - M(:, 2)];
SCALEFACTOR = 3;   % use 3-sigma rule

[srcLength temp0] = size(srcPoints);
[tgtLength temp0] = size(tgtPoints);

% Step 1: Compute the distance matrix D. 
% Notice that this can also be done by creating a C-Mex function similar to createDistanceMatrix

D = zeros(tgtLength, srcLength, 2); 
% Set dimension to 2 to allow computation of distance along and perpendicular to flow vector

for i = 1 : tgtLength
    flowVecMag = sqrt(flowVecList(i, 1)^2 + flowVecList(i, 2)^2);
    % This step can be omitted if the vectors have been normalized.
    for j = 1 : srcLength
        relativeD = [tgtPoints(i, 1) - srcPoints(j, 1), tgtPoints(i, 2) - srcPoints(j, 2)];
        relativeDMag = sqrt(relativeD(1)^2 + relativeD(2)^2);  
        
        if (relativeDMag < 1.e-6)
            D(i, j, 1) = 0;
            D(j, i, 2) = 0;  % Must prevent dividing by zero
        else
            theta = acos(relativeD * flowVecList(i, :)' / relativeDMag / flowVecMag);
            % This angle is actually coordinate-system-independent.
            
            temp = relativeDMag * [cos(theta) sin(theta)];
            if (temp(1) < SCALEFACTOR * sigmaU) 
                D(i, j, 1) = temp(1);
            else
                D(i, j, 1) = -1; % -1 indicates that the point is not within the search region
            end
        
            if (temp(2) < SCALEFACTOR * sigmaV) 
                D(i, j, 2) = temp(2);
            else
                D(i, j, 2) = -1;
            end
        end
        
        if (i ~= j)
            D(j, i, :) = D(i, j, :);
        end
    end
end

% Step 2: Computer the flow vector for each point
for i = 1 : tgtLength
    % Generate the list of points within the given radius
    weightsum = 0;
    vecsum = [0 0];
    % Calculate the weight
    for j = 1 : srcLength
        if ((D(i, j, 1) ~= -1) & (D(i, j, 2) ~= -1))
            weight = exp(-0.5 * D(i, j, 1)^2 / (sigmaU^2) - 0.5 * D(i, j, 2)^2 / (sigmaV^2));
            % since normalization will be done, it is not necessary to multiply using 0.5 / pi / sigmaU / sigmaV *
            weightsum = weightsum + weight;
            vecsum = vecsum + weight * V(j, :);
        end
    end
    interField(i, 1 : 2) = vecsum / weightsum;  % Normalization
end





