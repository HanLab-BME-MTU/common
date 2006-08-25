function [numObsPerBin,binCenter,gaussParam,errFlag] = fitHistWithGaussians(...
    observations,alpha,variableMean,variableVar,showPlot)
%FITHISTWITHGAUSSIANS fits multiple Gaussians to a histogram (including determining the number of necessary Gaussians)
%
%SYNOPSIS [numObsPerBin,binCenter,gaussParam,errFlag] = fitHistWithGaussians(...
%    observations,alpha,variableMean,variableVar,showPlot)
%
%INPUT  observations: Vector of observations whose histogram is to be fitted.
%       alpha       : Alpha-value for the statistical test that compares the 
%                     fit of n+1 Gaussians to the fit of n Gaussians. 
%                     Optional. Default: 0.05.
%       variableMean: 0 if assuming the fixed relationship
%                     (mean of nth Gaussian) = n * (mean of 1st Gaussian).
%                     1 if there is no relationship between the means of 
%                     different Gaussians. 
%                     Optional. Default: 0.
%       variableVar : 0 if assuming that all Gaussians have the same
%                     variance. 1 if there is no relationshop between the 
%                     variances of different Gaussians. 
%                     Optional. Default: 0.
%       showPlot    : 1 to plot the histogram and fitted Gaussians, 0 otherwise. 
%                     Optional. Default: 1.
%
%OUTPUT numObsPerBin: Number of observations that fall in each bin.
%       binCenter   : Center of each bin.
%       gaussParam  : Matrix with number of rows equal to number of fitted
%                     Gaussians and three columns indicating the mean,
%                     variance and amplitude of each Gaussian.
%       errFlag     : 0 if function executes normally, 1 otherwise.
%
%REMARKS The fitted Gaussians are normalized. Thus, the contribution of one
%Gaussian is given by
%(gaussParam(3)/(gaussParam(2)*sqrt(2pi)))
%                     *exp(-(x-gaussParam(1))^2/(2*gaussParam(2)^2)
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numObsPerBin = [];
binCenter = [];
gaussParam = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--fitHistWithGaussians: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(alpha)
    alpha = 0.05;
else
    if alpha <= 0 || alpha >= 1
        disp('fitHistWithGaussians: Variable "alpha" should be between 0 and 1!');
    end
end

if nargin < 3 || isempty(variableMean)
    variableMean = 0;
else
    if variableMean ~= 0 && variableMean ~= 1
        disp('fitHistWithGaussians: Variable "variableMean" should be either 0 and 1!');
    end
end

if nargin < 4 || isempty(variableVar)
    variableVar = 0;
else
    if variableVar ~= 0 && variableVar ~= 1
        disp('fitHistWithGaussians: Variable "variableVar" should be either 0 and 1!');
    end
end

if nargin < 5 || isempty(showPlot)
    showPlot = 1;
else
    if showPlot ~= 0 && showPlot ~= 1
        disp('fitHistWithGaussians: Variable "showPlot" should be either 0 and 1!');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Histogram calculation and fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the number of observations
numObservations = length(find(~isnan(observations)));

%calculate the histogram
[numObsPerBin,binCenter] = histogram(observations);
numObsPerBin = numObsPerBin';
binCenter = binCenter';

%determine the number of bins used
numBins = length(binCenter);

%calculate the cumulative histogram
cumHist = zeros(numBins,1);
for iBin = 1 : numBins
    cumHist(iBin) = sum(numObsPerBin(1:iBin));
end

%initialize variables indicating number of fitted Gaussians and their parameters
numGauss = 0;
gaussParam = [];

%logical variable indicating whether to attempt to fit
fit = 1;

%set some optimization options
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-6);

%fit the cumulative histogram with as many Gaussians as necessary
while fit

    if variableMean %if mean is variable
        
        if variableVar %if variance is variable
            
            %add another Gaussian to the fit
            numGaussT = numGauss + 1;

            %calculate number of degrees of freedom
            numDegFreeT = numBins - 3*numGaussT;

            %assign initial values to unknown parameters
            gaussParamT = [gaussParam; [binCenter(floor(end/2)) ...
                10*(binCenter(end)-binCenter(end-1)) numObservations]];
            x0 = gaussParamT(:);

            %assign lower bounds
            lb = [binCenter(1)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(end-1))*ones(numGaussT,1) ...
                zeros(numGaussT,1)];
            lb = lb(:);
            
            %assign upper bounds
            ub = [binCenter(end)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(1))*ones(numGaussT,1) ...
                numObservations*ones(numGaussT,1)];
            ub = ub(:);
            
            %estimate unknown parameters
            [param,resnorm,residualsT] = lsqcurvefit(@calcCumDistrNGauss,x0,...
                binCenter,cumHist,lb,ub,options,variableMean,variableVar);

            %get output from parameters vector
            gaussParamT = reshape(param,numGaussT,3);
            
        else %if variance is constrained

            %add another Gaussian to the fit
            numGaussT = numGauss + 1;

            %calculate number of degrees of freedom
            numDegFreeT = numBins - 2*numGaussT - 1;

            %assign initial values to unknown parameters
            gaussParamT = [gaussParam; [binCenter(floor(end/2)) ...
                10*(binCenter(end)-binCenter(end-1)) numObservations/2]];
            x0 = [gaussParamT(:,1); gaussParamT(1,2); gaussParamT(:,3)];

            %assign lower bounds
            lb = [binCenter(1)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(end-1))*ones(numGaussT,1) ...
                zeros(numGaussT,1)];
            lb = [lb(:,1); lb(1,2); lb(:,3)];

            %assign upper bounds
            ub = [binCenter(end)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(1))*ones(numGaussT,1) ...
                numObservations*ones(numGaussT,1)];
            ub = [ub(:,1); ub(1,2); ub(:,3)];

            %estimate unknown parameters
            [param,resnorm,residualsT] = lsqcurvefit(@calcCumDistrNGauss,x0,...
                binCenter,cumHist,lb,ub,options,variableMean,variableVar);

            %get output from parameters vector
            gaussParamT(:,1) = param(1:numGaussT);
            gaussParamT(:,2) = repmat(param(numGaussT+1),numGaussT,1);
            gaussParamT(:,3) = param(numGaussT+2:end);

        end %(if variableVar ... else ...)

    else %if mean is constrained

        if variableVar %if variance is variable
            
            %add another Gaussian to the fit
            numGaussT = numGauss + 1;

            %calculate number of degrees of freedom
            numDegFreeT = numBins - 2*numGaussT - 1;

            %assign initial values to unknown parameters
            gaussParamT = [gaussParam; [binCenter(floor(end/2)) ...
                10*(binCenter(end)-binCenter(end-1)) numObservations/2]];
            x0 = [gaussParamT(1,1); gaussParamT(:,2); gaussParamT(:,3)];

            %assign lower bounds
            lb = [binCenter(1)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(end-1))*ones(numGaussT,1) ...
                zeros(numGaussT,1)];
            lb = [lb(1,1); lb(:,2); lb(:,3)];

            %assign upper bounds
            ub = [binCenter(end)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(1))*ones(numGaussT,1) ...
                numObservations*ones(numGaussT,1)];
            ub = [ub(1,1); ub(:,2); ub(:,3)];

            %estimate unknown parameters
            [param,resnorm,residualsT] = lsqcurvefit(@calcCumDistrNGauss,x0,...
                binCenter,cumHist,lb,ub,options,variableMean,variableVar);

            %get output from parameters vector
            gaussParamT(:,1) = [1:numGaussT]'*param(1);
            gaussParamT(:,2) = param(2:numGaussT+1);
            gaussParamT(:,3) = param(numGaussT+2:end);

        else %if variance is constrained

            %add another Gaussian to the fit
            numGaussT = numGauss + 1;

            %calculate number of degrees of freedom
            numDegFreeT = numBins - numGaussT - 2;

            %assign initial values to unknown parameters
            gaussParamT = [gaussParam; [binCenter(floor(end/2)) ...
                10*(binCenter(end)-binCenter(end-1)) numObservations/2]];
            x0 = [gaussParamT(1,1:2)'; gaussParamT(:,3)];

            %assign lower bounds
            lb = [binCenter(1)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(end-1))*ones(numGaussT,1) ...
                zeros(numGaussT,1)];
            lb = [lb(1,1:2)'; lb(:,3)];

            %assign upper bounds
            ub = [binCenter(end)*ones(numGaussT,1) ...
                (binCenter(end)-binCenter(1))*ones(numGaussT,1) ...
                numObservations*ones(numGaussT,1)];
            ub = [ub(1,1:2)'; ub(:,3)];

            %estimate unknown parameters
            [param,resnorm,residualsT] = lsqcurvefit(@calcCumDistrNGauss,x0,...
                binCenter,cumHist,lb,ub,options,variableMean,variableVar);

            %get output from parameters vector
            gaussParamT(:,1) = [1:numGaussT]'*param(1);
            gaussParamT(:,2) = repmat(param(2),numGaussT,1);
            gaussParamT(:,3) = param(3:end);

        end %(if variableVar ... else ...)

    end %(if variableMean ... else ...)

    %check whether addition of 1 Gaussian has significantly improved the fit
    if numGaussT > 1 %if this is not the first fit

        %get test statistic, which is F-distributed
        testStat = (sum(residualsT.^2)/numDegFreeT)/...
            (sum(residuals.^2)/numDegFree);

        %get p-value of test statistic
        pValue = fcdf(testStat,numDegFreeT,numDegFree);

        %compare p-value to alpha
        %1-sided F-test: H0: F=1, H1: F<1
        if pValue < alpha %if p-value is smaller
            fit = 1; %accept this fit and attempt another fit with an additional Gaussian
        else %if p-value is larger
            fit = 0; %do not accept this fit and exit
        end

    end %(if numGaussT > 1)

    %if this fit is accepted, update some variables
    if fit
        numGauss = numGaussT;
        gaussParam = gaussParamT;
        residuals = residualsT;
        numDegFree = numDegFreeT;
    end

end %(while fit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if the user wants to plot
if showPlot

    %get the distribution from the optimized parameters
    distrNGauss = zeros(size(binCenter));
    for i=1:numGauss
        distrNGauss = distrNGauss + gaussParam(i,3)*normpdf(binCenter,...
            gaussParam(i,1),gaussParam(i,2))*(binCenter(end)-binCenter(end-1));
    end

    %get the cumulative distribution from the optimized parameters
    cumDistrNGauss = zeros(size(binCenter));
    for i=1:numGauss
        cumDistrNGauss = cumDistrNGauss + gaussParam(i,3)*normcdf(binCenter,...
            gaussParam(i,1),gaussParam(i,2));
    end

    %make new figure
    figure

    %plot the histogram and the fitted Gaussians in the left half of the
    %figure
    subplot(1,2,1);
    plot(binCenter,numObsPerBin,'k.')
    hold on
    plot(binCenter,distrNGauss,'r')

    %plot the histogram and the fitted Gaussians in the right half of the
    %figure
    subplot(1,2,2);
    plot(binCenter,cumHist,'k.')
    hold on
    plot(binCenter,cumDistrNGauss,'r')

end %(if showPlot)

%%%%% ~~ the end ~~ %%%%%
