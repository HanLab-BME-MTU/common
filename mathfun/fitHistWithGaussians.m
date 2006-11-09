function [numObsPerBinP,binCenterP,gaussParam,errFlag] = fitHistWithGaussians(...
    observations,alpha,variableMean,variableVar,showPlot,maxNumGauss)
%FITHISTWITHGAUSSIANS fits multiple Gaussians to a histogram (including determining the number of necessary Gaussians)
%
%SYNOPSIS [numObsPerBin,binCenter,gaussParam,errFlag] = fitHistWithGaussians(...
%    observations,alpha,variableMean,variableVar,showPlot,maxNumGauss)
%
% or [numObsPerBin,binCenter,gaussParam,errFlag] = fitHistWithGaussians(...
%    observations,'R',showPlot,maxNumGauss)
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
%       showPlot    : 0 to not plot anything
%                     1 to plot the histogram and fitted Gaussians
%                     2 as 1, but with smooth histogram
%                     Optional. Default: 1.
%       maxNumGauss : upper limit to the number of Gaussians.
%                         Optional. Default: 100 (if 'R', default: 9)
%       'R'         : If you input 'R' as second element, the program will
%                     call the function mclust in the statistical package
%                     R. For this, you will need the toolbox matlab2R, and
%                     a local installation of R with a COM interface, and
%                     Windows as OS.
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
%       The function does not yet allow to set min#clust, maxPoints!
%
%Khuloud Jaqaman, August 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numObsPerBinP = [];
binCenterP = [];
gaussParam = [];
errFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1 || isempty(observations)
    disp('--fitHistWithGaussians: Incorrect number of input arguments!');
    errFlag = 1;
else
    % make sure observations is a col-vector
    observations = observations(:);
end

% Check whether to use the matlab routine or R
switch isnumeric(alpha)
    case 1

        isR = false;

        if nargin < 2 || isempty(alpha)
            alpha = 0.05;
        else
            if alpha <= 0 || alpha >= 1
                disp('--fitHistWithGaussians: Variable "alpha" should be between 0 and 1!');
                errFlag = 1;
            end
        end

        if nargin < 3 || isempty(variableMean)
            variableMean = 0;
        else
            if variableMean ~= 0 && variableMean ~= 1
                disp('--fitHistWithGaussians: Variable "variableMean" should be either 0 and 1!');
                errFlag = 1;
            end
        end

        if nargin < 4 || isempty(variableVar)
            variableVar = 0;
        else
            if variableVar ~= 0 && variableVar ~= 1
                disp('--fitHistWithGaussians: Variable "variableVar" should be either 0 and 1!');
                errFlag = 1;
            end
        end

        if nargin < 5 || isempty(showPlot)
            showPlot = 1;
        else
            if ~any(showPlot == [0,1,2])
                disp('--fitHistWithGaussians: Variable "showPlot" should be either 0, 1, or 2!');
                errFlag = 1;
            end
        end

        if nargin < 6 || isempty(maxNumGauss)
            maxNumGauss = 100;
        else
            if maxNumGauss < 1
                disp('--fitHistWithGaussians: Variable "maxNumGaussians" should be at least 1');
                errFlag = 1;
            end
        end

    case 0 % alpha is not numeric

        isR = true;

        if strmatch(alpha,'R')
            % check whether we are on windows, and try to launch R
            if ispc
                try
                    status = openR;
                catch
                    disp('--fitHistWithGaussians: unable to launch R');
                    errFlag = 1;
                end
                if ~status
                    if strcmp(lastwarn,'Already connected to an R server.')
                        leaveRopen = true;
                    else
                        disp('--fitHistWithGaussians: unable to launch R');
                        errFlag = 1;
                    end
                else
                    leaveRopen = false;
                end
            else
                disp('--fitHistWithGaussians: ''R'' needs the COM interface and therefore only runs under Windows');
                errFlag = 1;
            end
        else
            disp('--fitHistWithGaussians: unknown option for second input');
            errFlag = 1;
        end

        % check for display and maxNumGaussians
        % variableMean ~showPlot
        if nargin < 3 || isempty(variableMean)
            showPlot = 1;
        else
            if ~any(variableMean == [0,1,2])
                disp('--fitHistWithGaussians: Variable "showPlot" should be either 0, 1, or 2!');
                errFlag = 1;
            else
                showPlot = variableMean;
            end
        end
        % variableVar ~maxNumGauss
        if nargin < 4 || isempty(variableVar)
            maxNumGauss = 9;
        else
            if variableVar < 1
                disp('--fitHistWithGaussians: Variable "maxNumGaussians" should be at least 1');
                errFlag = 1;
            else
                maxNumGauss = variableVar;
            end
        end

end % switch

% check error flag
if errFlag
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram calculation and fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code bit replaced to avoid any problems with binning
% %get the number of observations
% numObservations = length(find(~isnan(observations)));
%
% %calculate the histogram
% [numObsPerBin,binCenter] = histogram(observations);
% numObsPerBin = numObsPerBin';
% binCenter = binCenter';
%
% %determine the number of bins used
% numBins = length(binCenter);
%
% %calculate the cumulative histogram
% cumHist = zeros(numBins,1);
% for iBin = 1 : numBins
%     cumHist(iBin) = sum(numObsPerBin(1:iBin));
% end


% ---------- SWITCH BETWEEN MATLAB AND R ---------------

switch isR

    case 0 % run Matlab

        % for the optimization: don't bin the cumulative histogram. However, don't
        % use duplicate values - therefore, use cdfcalc. It also returns the number
        % of non-NaN observations, and an error message, if any.
        [cumHist,binCenter,numObservations,errMsg] = cdfcalc(observations);
        if ~isempty(errMsg)
            % disp/return instead of throwing the error b/c of Khuloud's standard
            disp(sprintf('--%s',errMsg))
            return
        end
        % number of bins is the number of different x-values
        numBins = length(binCenter);
        % cdfcalc returns n+1 values for cumHist. 1:end-1 is the bottom of the
        % step, 2:end the top. Take the middle for best results.
        % cumHist = (cumHist(2:end)+cumHist(1:end-1))/2;

        % make cumHist with binCenters in middle of top of step
        binCenter = (binCenter(1:end-1)+binCenter(2:end))/2;
        cumHist = cumHist(2:end-1);
        numBins = numBins - 1;


        % downsample to about 1000 points if necessary
        if numBins > 1000
            dsIdx = unique(round(linspace(1,numBins,1000)))';
            cumHist = cumHist(dsIdx);
            binCenter = binCenter(dsIdx);
            numBins = length(dsIdx);
        end

        % make cumHist go from 1:numObservations
        cumHist = cumHist * numObservations;

        %initialize variables indicating number of fitted Gaussians and their parameters
        numGauss = 0;
        gaussParam = [];

        %logical variable indicating whether to attempt to fit
        fit = 1;

        %set some optimization options
        options = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-6,'Display','off');

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
                    gaussParamT(:,1) = (1:numGaussT)'*param(1);
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
                    gaussParamT(:,1) = (1:numGaussT)'*param(1);
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
                if pValue < alpha && numGauss <= maxNumGauss %if p-value is smaller and the limit of Gaussians isn't reached
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

        % ----- R ------
    case 1
        % check if downsampling necessary - cap at 1000 observations
        numObservations = length(observations);
        % remove NaN, Inf
        observationsC = observations;
        observationsC(~isfinite(observationsC)) = [];
        numObservationsC = length(observationsC);
        if  numObservationsC > 1000
            observationsDS = sort(observationsC);
            observationsDS = observationsDS(round(linspace(1,numObservationsC,1000)));
        else
            observationsDS = observationsC;
        end
        clear observationsC

        % run R (it's already open)

        % load mclust
        evalR('library(mclust)')

        % put observations into R
        putRdata('inputData',observationsDS);

        % run Mclust - they have changed the syntax since the last version!
        try
            % new version
        evalR(sprintf('clusterData <- Mclust(inputData,G=1:%i)',maxNumGauss))
        % read mu, sigma, and probabilities (which will give amplitudes
        % when multiplied by the number of observations)
        mu = evalR('clusterData$parameters$mean');
        sigma = sqrt(evalR('clusterData$parameters$variance$sigmasq'));
        numGauss = length(mu);
        if numGauss == 1
            amp = numObservations;
        else
        amp = evalR('clusterData$parameters$pro')*numObservations;
        end
        
        catch
            try
                % old version
                evalR(sprintf('clusterData <- Mclust(inputData,1,%i)',maxNumGauss))
                % read mu, sigma, and probabilities (which will give amplitudes
        % when multiplied by the number of observations)
        mu = evalR('clusterData$mu');
        sigma = sqrt(evalR('clusterData$sigma'));
        % take care of possible equal sigma
        if length(sigma) == 1
            sigma = repmat(sigma,1,length(mu));
        end
        amp = evalR('clusterData$pro')*numObservations;
        numGauss = length(mu);
            catch
                rethrow(lasterror)
            end
        end

        % check the number of sigmas
        if length(sigma) ~= numGauss
            sigma = sigma * ones(1,numGauss);
        end

        % catenate into gaussParam
        gaussParam = [mu',sigma',amp'];



        % close R
        if ~leaveRopen
            closeR
        end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make a histogram for plotting and output. Choose how to calculate bins
if showPlot == 2
    [numObsPerBinP,binCenterP] = histogram(observations,'smooth');
else
    [numObsPerBinP,binCenterP] = histogram(observations);
end

%if the user wants to plot
if showPlot

    %get the distribution from the optimized parameters
    distrNGauss = zeros(size(binCenterP));
    for i=1:numGauss
        % no longer multiply by the bin width - histograms is normed now
        distrNGauss = distrNGauss + gaussParam(i,3)*normpdf(binCenterP,...
            gaussParam(i,1),gaussParam(i,2)) ;%*(binCenterP(end)-binCenterP(end-1));
    end

    if isR
        [cumHist,binCenter] = cdfcalc(observations);
        binCenter = (binCenter(1:end-1)+binCenter(2:end))/2;
        cumHist = cumHist(2:end-1) * numObservations;
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
    %figure. Correct by the number of NaNs
    subplot(1,2,1);
    bar(binCenterP,numObsPerBinP,'k')
    hold on
    plot(binCenterP,distrNGauss * sum(isfinite(observations))/numObservations,'r')

    %plot the histogram and the fitted Gaussians in the right half of the
    %figure
    subplot(1,2,2);
    plot(binCenter,cumHist,'k.')
    hold on
    plot(binCenter,cumDistrNGauss,'r')

end %(if showPlot)

%%%%% ~~ the end ~~ %%%%%
