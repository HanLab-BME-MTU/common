function asymParam = asymDeterm2D3D(positions)
%ASYMDETERM2D3D estimates the asymmetry in a scatter of positions
%
%SYNPOSIS asymParam = asymDeterm2D3D(positions,probDim)
%
%INPUT  positions: n-by-2/3 array of positions (x,y,[z])
%
%OUTPUT asymParam: Parameter estimating asymmetry of positional scatter
%
%Khuloud Jaqaman, December 2007

%get problem dimensionality (2D or 3D)
probDim = size(positions,2);

%calculate the variance-covariance matrix of positions
posCov = nancov(positions);

%perform eigenvalue decomposition
[eigenVec,eigenVal] = eig(posCov);

%get the eigenvalues in a vector
eigenVal = diag(eigenVal);

%calculate some intermediate sums
doubleSum = 0;
for i = 1 : probDim - 1
    doubleSum = doubleSum + sum( ( eigenVal(i) - eigenVal(i+1:end) ).^2 );
end
singleSum = (probDim - 1) * ( sum(eigenVal) )^2;

%calculate asymmetry parameter
asymParam = -log( 1 - doubleSum / singleSum );
