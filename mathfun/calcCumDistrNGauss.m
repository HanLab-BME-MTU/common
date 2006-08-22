function cumDistrNGauss = calcCumDistrNGauss(param,abscissa,variableMean,...
    variableVar)
%CALCCUMDISTRNGAUSS calculates the cumulative distribution of N Gaussians
%
%SYNOPSIS cumDistrNGauss = calcCumDistrNGauss(param,abscissa,variableMean,...
%    variableVar)
%
%INPUT  param         : Vector of parameters indicating the means,
%                       variances and amplitudes of the N Gaussians.
%                       -If variableMean=1 & variableVar=1, param has 3N 
%                        entries: N means, N variances and N amplitudes.
%                       -If variableMean=1 & variableVar=0, param has 2N+1
%                        entries: N means, 1 variance and N amplitudes.
%                       -If variableMean=0 & variableVar=1, param has 2N+1
%                        entries: 1 mean, N variances and N amplitudes.
%                       -If variableMean=0 & variableVar=0, param has N+2
%                        entries: 1 mean, 1 variance and N amplitudes.
%                       See below the definitions of variableMean and
%                       variableVar.
%       abscissa      : Abscissa values at which the cumulative
%                       distribution is calculated.
%       variableMean  : 0 if assuming the fixed relationship
%                       (mean of nth Gaussian) = n * (mean of 1st Gaussian).
%                       1 if there is no relationship between the means of 
%                       different Gaussians. 
%                       Optional. Default: 1.
%       variableVar   : 0 if assuming that all Gaussians have the same
%                       variance. 1 if there is no relationshop between the 
%                       variances of different Gaussians. 
%                       Optional. Default: 1.
%
%OUTPUT cumDistrNGauss: Values of the resulting cumulative distribution
%                       given the input abscissa values.

%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cumDistrNGauss = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 2
    disp('--calcCumDistrNGauss: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

if nargin < 3 || isempty(variableMean)
    variableMean = 1;
end

if nargin < 4 || isempty(variableVar)
    variableVar = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the cumulative distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the means, variances and amplitudes of the Gaussians from the input
%parameter vector
if variableMean %if mean is variable

    if variableVar %if variance is variable

        %get number of Gaussians
        numGauss = length(param)/3;

        %get their means, variances and amplitudes
        gaussMean = param(1:numGauss);
        gaussVar  = param(numGauss+1:2*numGauss);
        gaussAmp  = param(2*numGauss+1:end);

    else %if variance is not variable

        %get number of Gaussians
        numGauss = floor(length(param)/2);

        %get their means, variances and amplitudes
        gaussMean = param(1:numGauss);
        gaussVar  = repmat(param(numGauss+1),numGauss,1);
        gaussAmp  = param(numGauss+2:end);

    end %(if variableVar ... else ...)

else %if mean is not variable

    if variableVar %if variance is variable

        %get number of Gaussians
        numGauss = floor(length(param)/2);

        %get their means, variances and amplitudes
        gaussMean = [1:numGauss]'*param(1);
        gaussVar  = param(2:numGauss+1);
        gaussAmp  = param(numGauss+2:end);

    else %if variance is not variable

        %get number of Gaussians
        numGauss = length(param)-2;

        %get their means, variances and amplitudes
        gaussMean = [1:numGauss]'*param(1);
        gaussVar  = repmat(param(2),numGauss,1);
        gaussAmp  = param(3:end);

    end %(if variableVar ... else ...)

end %(if variableMean ... else ...)

%calculate the cumulative distribution
cumDistrNGauss = zeros(size(abscissa));
for i=1:numGauss
    cumDistrNGauss = cumDistrNGauss + gaussAmp(i)*normcdf(abscissa,...
        gaussMean(i),gaussVar(i));
end

%%%%% ~~ the end ~~ %%%%%

