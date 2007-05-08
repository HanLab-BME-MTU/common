function asymmetry = asymDetermination(positions);

% ASYMDETERMINATION  function which calculates the asymmetry of the gyration
%                    tensor of the x and y values of a given set of
%                    datapoints (trajectory)
%
% SYNOPSIS  asymmetry = asymDetermination(positions)
% 
% INPUT     positions = n-by-m matrix
%                   
%                   rows: timepoints along a trajectory
%                   columns: 1. x-values
%                            2. y-values  
%
% OUTPUT    asymmetry = asymmetry of all positions in the dataset
%                       (trajectory)
%
% CREATED gp 4/02/07

%-----------------------------------------
% Initialize
%-----------------------------------------

asymmetry = [];

%----------------------------------------------------------
% Calculation of the 2D-radius (Rg) of the gyration tensor 
%----------------------------------------------------------

R11 = max((mean(positions(:,1).^2)) - ((mean(positions(:,1))).^2),eps);

R22 = max((mean(positions(:,2).^2)) - ((mean(positions(:,2))).^2),eps);

R12 = max((mean(positions(:,1).*positions(:,2))) - ...
    (mean(positions(:,1)) * (mean(positions(:,2)))),eps);

R21 = R12;

Rg = [R11 R12; R21 R22];

% computes the eigen value and the eigenvector of gyration tensor Rg
[eigenVector, eigenValue] = eig(Rg);

Rxx = sqrt(eigenValue(1,1));

Ryy = sqrt(eigenValue(2,2));

% calculation of the asymmetry 
asymmetry = -log(1-(((Rxx^2 - Ryy^2)^2) / (Rxx^2 + Ryy^2)^2));