function [objFunc,residuals] = multiStepFunction(param,x,y)

%Khuloud Jaqaman, July 2014

%% Output

objFunc = [];
residuals = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--multiStepFunction: Incorrect number of input arguments!');
    return
end

if nargin < 3 || isempty(y)
    y = [];
end

%% Calculation

%get number of steps
numSteps = (length(param)-1)/2;

%get step locations
stepPos = param(1:numSteps);

%get function values between steps
betweenStepsVal = param(numSteps+1:end);

%go over x and calculate function
multiStepVal = betweenStepsVal(1)*ones(size(x));
for iStep = 1 : numSteps
    multiStepVal(x>stepPos(iStep)) = betweenStepsVal(iStep+1);
end

%take difference with data and calculate objective function
if isempty(y)
    residuals = [];
    objFunc = [];
else
    residuals = multiStepVal - y;
    objFunc = residuals' * residuals;
end


%% ~~~ the end ~~~

