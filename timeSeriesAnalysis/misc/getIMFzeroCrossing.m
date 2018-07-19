function zC = getIMFzeroCrossing(imf)
%This function calculates the number of zero-crossing of the input
%
%Input:
%       imf - #ofVariables,#ofPoints
%
%Marco Vilela
%2013

[~,nPoint] = size(imf);

%Number of zero-crossing
zC = sum( imf(:,1:nPoint-1).*imf(:,2:nPoint) < 0,2);

end

