function asymmetry = asymDeterm2D3D(positions,probDim)

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
%                            3. z-values (if 3D)
%
%           probDim: Problem dimension. 2 (for 2D) or 3 (for 3D)
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

%reserve memory for the gyration tensor Rg
Rg = zeros(probDim);

%calculate the components of Rg
for i = 1 : probDim
    for j = i : probDim
        Rg(i,j) = max((mean(positions(:,i).*positions(:,j))) - ...
            (mean(positions(:,i)) * (mean(positions(:,j)))),eps);
        Rg(j,i) = Rg(i,j); %tensor is symmetric
    end
end

%compute the eigenvalues and the eigenvectors of Rg
[eigenVector,eigenValue] = eig(Rg);

%get the eigenvalues in a vector
eigenValue = diag(eigenValue);

%calculate the asymmetry parameter
doubleSum = 0;
for i = 1 : probDim - 1
    doubleSum = doubleSum + sum((eigenValue(i) - eigenValue(i+1:end)).^2);
end
singleSum = (probDim - 1) * (sum(eigenValue))^2;
asymmetry = -log(1 - doubleSum / singleSum);

