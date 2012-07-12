function [diffModeDivider,fracTruePos] = estimDiffModeDividers(diffModeCoef,numTraj,trajLength,doPlot)
%ESTIMDIFFMODEDIVIDERS estimates the optimal dividers between different diffusion modes
%
%SYNOPSIS [diffModeDivider,fracTruePos] = estimDiffModeDividers(diffModeCoef,numTraj,trajLength,doPlot)
%
%INPUT  diffModeCoef: A row vector listing the diffusion coefficient of the
%                     different diffusion modes. The diffusion modes should
%                     be sorted in ascending order of their diffusion
%                     coefficients.
%       numTraj     : Number of trajectories to be simulated.
%                     Optional. Default: 1000.
%       trajLength  : Number of time points per trajectory.
%                     Optional. Default: 20.
%       doPlot      : 1 to plot the results, 0 otherwise.
%                     Optional. Default: 0.
%
%OUTPUT diffModeDivider: Estimated dividers between the different modes.
%       fracTruePos    : Fraction of each mode correctly classified using
%                        the estimated dividers.
%
%REMARKS Code is written for 2D case only, but can be generalized to 3D.
%
%Khuloud Jaqaman, July 2012

%% Output
diffModeDivider = [];
fracTruePos = [];

%% Input

if nargin < 1
    disp('estimDiffModeDividers: Please enter vector of diffusion modes')
    return
end

if nargin < 2 || isempty(numTraj)
    numTraj = 1000;
end

if nargin < 3 || isempty(trajLength)
    trajLength = 20;
end

if nargin < 4 || isempty(doPlot)
    doPlot = 0;
end

%% Divider estimation

%get number of diffusion modes
numMode = length(diffModeCoef);

%simulate trajectories in each mode
traj = NaN(trajLength,2,numTraj,numMode);
for iMode = 1 : numMode
    for iTraj = 1 : numTraj
        trajTmp = brownianMotion(2,diffModeCoef(iMode),trajLength,0.01);
        trajTmp = trajTmp(1:100:end,:);
        trajTmp = trajTmp(2:end,:);
        traj(:,:,iTraj,iMode) = trajTmp;
    end
end
msd = squeeze( mean( sum( diff(traj).^2,2 ) ) );
diffCoef = msd / 4;

%for each pair of consecutive modes, slide divider between their means and
%calculate the fraction of each mode that will be properly classified
divValue = 0:0.001:diffModeCoef(end)*1.1;
numValue = length(divValue);
fracModeBelow = NaN(numValue,numMode);
fracModeAbove = NaN(numValue,numMode);
for iMode = 1 : numMode
    for iValue = 1 : numValue
        fracModeBelow(iValue,iMode) = length(find(diffCoef(:,iMode)<=divValue(iValue)))/numTraj;
        fracModeAbove(iValue,iMode) = length(find(diffCoef(:,iMode)>divValue(iValue)))/numTraj;
    end
end

%find the best divider between each pair of modes, i.e. the divider value
%which simultaneously maximizes the correct classification of each mode
diffModeDivider = NaN(numMode-1,1);
fracTruePos = NaN(numMode-1,1);
for iMode = 1 : numMode-1
    [dividerTmp,fracTmp] = polyxpoly(divValue',fracModeBelow(:,iMode),divValue',fracModeAbove(:,iMode+1));
    diffModeDivider(iMode) = mean(dividerTmp);
    fracTruePos(iMode) = fracTmp(1);
end

%for display, remove unnecessary regions
fracModeBelowDisp = fracModeBelow;
fracModeAboveDisp = fracModeAbove;
for iMode = 1 : numMode
    fracModeBelowDisp(divValue<diffModeCoef(iMode),iMode) = NaN;
    fracModeAboveDisp(divValue>diffModeCoef(iMode),iMode) = NaN;
    indx = find(fracModeBelowDisp(:,iMode)==1,6,'first');
    if ~isempty(indx)
        fracModeBelowDisp(indx(end):end,iMode) = NaN;
    end
    indx = find(fracModeAboveDisp(:,iMode)==1,6,'last');
    if ~isempty(indx)
        fracModeAboveDisp(1:indx(1),iMode) = NaN;
    end
end

%% Plotting

%plot the fractions and display dividers
if doPlot
    figure, hold on
    plot(divValue,fracModeBelowDisp(:,1:end-1))
    plot(divValue,fracModeAboveDisp(:,2:end))
    plot(diffModeDivider,fracTruePos,'ko')
end

%% ~~~ the end ~~~


% %display a histogram of the diffuion coefficients and indicate mode
% %dividers
% figure
% hist(log10(diffCoef),100)
% hold on
% for iMode = 1 : numMode-1
%     plot(log10(diffModeDivider(iMode))*[1 1],[0 250])
% end
%
% %calculate fraction of each mode included within dividers
% modeFracAboveBelowDivider = NaN(numMode,2);
% modeFracAboveBelowDivider(1,2) = length(find(diffCoef(:,1)<diffModeDivider(1)))/size(diffCoef,1);
% for i = 2 : numMode-1
%     modeFracAboveBelowDivider(i,1) = length(find(diffCoef(:,i)>diffModeDivider(i-1)))/size(diffCoef,1);
%     modeFracAboveBelowDivider(i,2) = length(find(diffCoef(:,i)<diffModeDivider(i)))/size(diffCoef,1);
% end
% i = numMode;
% modeFracAboveBelowDivider(i,1) = length(find(diffCoef(:,i)>diffModeDivider(i-1)))/size(diffCoef,1);
