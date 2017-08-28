function [polyParam,covMatBoot,polyParamBoot,xGrid,yOriginal,yBoot,yBootMean,yBootStd] = ...
    fitPolyBootstrapCI(dataX,dataY,polyOrder,bootSize,doPlot,stdsPlotCI)
%FITPOLYBOOTSTRAPCI fits a polynomial through input data points and bootstraps sample to esimate parameter variance-covariance and draw confidence interval around fit
%
%SYNPOSIS 
%   [polyParam,covMatBoot,polyParamBoot,xGrid,yOriginal,yBoot,yBootMean,yBootStd] = ...
%        fitPolyBootstrapCI(dataX,dataY,polyOrder,bootSize,doPlot,stdsPlotCI)
%
%INPUT  
%   Mandatory:
%       dataX, dataY: Column vectors of data to be fitted
%   Optional:
%       polyOrder   : Polynomial order. Default: 1.
%       bootSize    : Size of bootstrap sample Default: 1000.
%       doPlot      : 1 to plot data, fitted polynomial and confidence
%                     interval, 0 otherwise. Default: 1.
%       stdsPlotCI  : Number of stds to plot confidence interval. 
%                     1 x std gives 68% interval, 2 x std gives 95% interval,
%                     etc. Default: 1.
%
%OUTPUT
%       polyParam   : Row vector of fitted polynomial parameters, as output
%                     by polyfit.
%       covMatBoot  : Parameter variance-covariance matrix, as determined
%                     by bootstrapping. Order of parameters same as in
%                     polyParam.
%       polyParamBoot:Fitted polynomial parameters for bootstrap samples.
%                     Number of rows = number of bootstrap samples.
%       xGrid       : Row vector of x-values used for plotting fit. 5000
%                     values equally spaced from min(dataX) to max(dataX).
%       yOriginal   : Row vector of y-values from polynomial from fit to
%                     original data.
%       yBoot       : y-values from polynomial fit to boostrap samples.
%                     Number of rows = number of bootstrap samples.
%       yBootMean   : Row vector of mean y-value from all bootstrap
%                     samples.
%       yBootStd    : Row vector of y-value standard deviation from all
%                     bootstrap samples.
%
%Khuloud Jaqaman, June 2017

% %% Output
% 
% polyParam = [];
% covMatBoot = [];
% polyParamBoot = [];
% xGrid = [];
% yOriginal = [];
% yBoot = [];
% yBootMean = [];
% yBootStd = [];

%% Input

if nargin < 2
    error('Missing input data to be fitted!');
end

if nargin < 3 || isempty(polyOrder)
    polyOrder = 1;
end

if nargin < 4 || isempty(bootSize)
    bootSize = 1000;
end

if nargin < 5 || isempty(doPlot)
    doPlot = 1;
end

if nargin < 6 || isempty(stdsPlotCI)
    stdsPlotCI = 1;
end

%% Fitting and Bootstrapping

%fit polynomial to input data
polyParam = polyfit(dataX,dataY,polyOrder);

%bootstrap the polynomial fit
polyParamBoot = bootstrp(bootSize,@(x,y) polyfit(x,y,polyOrder),dataX,dataY);

%calculate parameter variance-covariance matrix
covMatBoot = cov(polyParamBoot);

%sample the x-axis values
xGrid = linspace(min(dataX),max(dataX),5000);

%make matrix of correct power to calculate y from x
xGridMat = ones(polyOrder+1,5000);
for iOrder = 1 : polyOrder
    xGridMat(polyOrder-iOrder+1,:) = xGrid.^iOrder;
end

%calculate y-values from x-values for fit to original data
yOriginal = polyParam * xGridMat;

%calculate y-values from x-values for each bootstrap sample
yBoot = zeros(bootSize,5000);
for iBoot = 1 : bootSize
    yBoot(iBoot,:) = polyParamBoot(iBoot,:) * xGridMat;
end

%calculate mean and std of y-values at each x-value
yBootMean = mean(yBoot);
yBootStd = std(yBoot);

%% Plotting

%plot if requested
if doPlot
    figure, hold on
    plot(dataX,dataY,'b.')
    plot(xGrid,yOriginal,'b')
    plot(xGrid,yBootMean,'b:')
    plot(xGrid,yBootMean+stdsPlotCI*yBootStd,'b--')
    plot(xGrid,yBootMean-stdsPlotCI*yBootStd,'b--')
end


