function [stepPos,betweenStepsVal] = fitDataWithMultiStepFunc(x,y,numStepsRange,alpha,doPlot)

%Khuloud Jaqaman, August 2014

%% Output
stepPos = [];
betweenStepsVal = [];

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

%get rid of NaNs in data
indx = find(~isnan(y));
x = x(indx);
y = y(indx);
numDataPoints = length(indx);

%% Calculation

%minimization options
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-3);

%fit models one by one, checking residuals
fit = 1;
firstFit = 1;
iStep = minNumSteps-1;
while fit && iStep < maxNumSteps
    
    %number of steps for current fit
    iStep = iStep + 1;
    
    %parameter initial guess
    stepPos0 = round(linspace(1,numDataPoints,iStep+2)');
    betweenStepsVal0 = ones(iStep+1,1);
    for iPart = 1 : iStep+1
        betweenStepsVal0(iPart) = mean(y(stepPos0(iPart):stepPos0(iPart+1)));
    end
    stepPos0 = x(stepPos0(2:end-1));
    
    %fit
    paramT = fminsearch(@multiStepFunction,[stepPos0; betweenStepsVal0],options,x,y);
    [~,residualsT] = multiStepFunction(paramT,x,y);
    numDegFreeT = numDataPoints - length(paramT);
    
    %keep fit if first, otherwise compare residuals
    if firstFit
        param = paramT;
        numDegFree = numDegFreeT;
        residuals = residualsT;
        firstFit = 0;
    else
        testStat = (sum(residualsT.^2)/numDegFreeT)/...
            (sum(residuals.^2)/numDegFree);
        pValue = fcdf(testStat,numDegFree,numDegFreeT);
        if pValue < alpha
            param = paramT;
            numDegFree = numDegFreeT;
            residuals = residualsT;
        else
            fit = 0;
        end
    end
    
end

%parameters for output
numSteps = (length(param)-1)/2;
stepPos = param(1:numSteps);
betweenStepsVal = param(numSteps+1:end);

%% Plotting

if doPlot
    
    figure, hold on
    plot(x,y)
    if numSteps == 0
        plot([x(1) x(end)],betweenStepsVal*[1 1],'r');
    else
        plot([x(1) stepPos(1)],betweenStepsVal(1)*[1 1],'r');
        plot(stepPos(1)*[1 1]',betweenStepsVal(1:2),'r');
        for iStep = 1 : numSteps-1
            plot([stepPos(iStep) stepPos(iStep+1)],betweenStepsVal(iStep+1)*[1 1],'r');
            plot(stepPos(iStep+1)*[1 1]',betweenStepsVal(iStep+1:iStep+2),'r');
        end
        plot([stepPos(end) x(end)],betweenStepsVal(end)*[1 1],'r');
    end
    
end


%% ~~~ the end ~~~

