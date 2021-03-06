function [stepX,valY] = fitDataWithMultiStepFunc2(x,y,numStepsRange,alpha,doPlot)
%FITDATAWITHMULTISTEPFUNC fits a time series with a multi-step function with number of steps determined automatically
%
%SYNOPSIS [stepX,valY] = fitDataWithMultiStepFunc(x,y,numStepsRange,alpha,doPlot)
%
%INPUT  x            : Independent variable of time series.
%       y            : Dependent variable of time series.
%       numStepsRange: Row vector indicating minimum and maximum number
%                      steps to attempt to fit.
%                      Optional. Default: [0 5].
%       alpha        : Alpha-value for F-test to determine number of steps.
%                      Optional. Default: 0.05.
%       doPlot       : 1 to plot data and fit, 0 otherwise.
%                      Optional. Default: 0.
%
%OUTPUT stepX        : x-values at which there is a step.
%       valY         : y-values between steps. If there are n steps, there
%                      will be n+1 y-values. First value is before first
%                      step, etc., and last value is after last step.
%
%Khuloud Jaqaman, August 2014

%% Output
stepX = [];
valY = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--fitDataWithMultiStepFunc: Incorrect number of input arguments!');
    return
end

if nargin < 3 || isempty(numStepsRange)
    numStepsRange = [0 5];
end

if nargin < 4 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 5 || isempty(doPlot)
    doPlot = 0;
end

%get range of steps to try to fit
minNumSteps = numStepsRange(1);
maxNumSteps = numStepsRange(2);

%take moving average of data to help with step initialization
[~,yMA] = movAvSmooth(y,5,2,0);

%get rid of NaNs in data
indx = find(~isnan(y));
x = x(indx);
y = y(indx);
yMA = yMA(indx);
numDataPoints = length(indx);

%% Calculation

%locate positions of largest jumps in data for step initialization
yMAdiff = diff(yMA);
[~,iSorted] = sort(yMAdiff);
stepPosGuess = x(iSorted);

%minimization options
options = optimset('MaxFunEvals',1000000,'MaxIter',100000,'TolFun',1e-4);

%fit models one by one, checking residuals
iTmp = 0;
ratioResid = NaN(maxNumSteps-minNumSteps+1,1);
for iStep = minNumSteps : maxNumSteps
    
    iTmp = iTmp + 1;
    
    %parameter initial guess
%     stepX0 = round(linspace(1,numDataPoints,iStep+2)');
    stepX0 = sort(stepPosGuess(1:iStep));
    stepIndx = [1; sort(iSorted(1:iStep)); numDataPoints];
    valY0 = ones(iStep+1,1);
    for iPart = 1 : iStep+1
        valY0(iPart) = mean(y(stepIndx(iPart):stepIndx(iPart+1)));
    end
    
    %fit
    paramT = fminsearch(@multiStepFunction,[stepX0; valY0],options,x,y);
    [~,tmp] = multiStepFunction(paramT,x,y);
    residualsT = tmp;
    stepX = paramT(1:iStep);
    valY = paramT(iStep+1:end);
    
    if iStep > 0
        
        %counter-fit
        stepX0 = [[x(1); stepX(1:end-1)] stepX];
        stepX0 = mean(stepX0,2);
        paramC = fminsearch(@multiStepFuncFS,valY,options,stepX0,x,y);
        [~,residualsC] = multiStepFuncFS(paramC,stepX0,x,y);
        
        %store ratio of residuals
        ratioResid(iTmp) = (residualsC'*residualsC)/(residualsT'*residualsT);
        
    end

end

%% Plotting

if doPlot
    
    overlayMultiStepFunc(x,y,stepX,valY)
    
%     figure, hold on
%     plot(x,y)
%     if numSteps == 0
%         plot([x(1) x(end)],valY*[1 1],'r');
%     else
%         plot([x(1) stepX(1)],valY(1)*[1 1],'r');
%         plot(stepX(1)*[1 1]',valY(1:2),'r');
%         for iStep = 1 : numSteps-1
%             plot([stepX(iStep) stepX(iStep+1)],valY(iStep+1)*[1 1],'r');
%             plot(stepX(iStep+1)*[1 1]',valY(iStep+1:iStep+2),'r');
%         end
%         plot([stepX(end) x(end)],valY(end)*[1 1],'r');
%     end
    
end


%% ~~~ the end ~~~

