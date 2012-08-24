function [modeParam,expParam,pvModeDiff,modeParamControl,pvModeDiffControl] = ...
    getDiffModes(tracksFinal,minLength,alpha,showPlot,maxNumMode,...
    binStrategy,plotName)
%GETDIFFMODES determines number of diffusion modes and their parameters from distribution of frame-to-frame displacements
%
%SYNOPSIS [modeParam,expParam,pvModeDiff,modeParamControl,pvModeDiffControl] = ...
%    getDiffModes(tracksFinal,minLength,alpha,showPlot,maxNumMode,...
%    binStrategy,plotName)
%
%INPUT  tracksFinal : Output of trackCloseGapsKalman.
%                     Optional. If not input, GUI will be called to get
%                     tracks. NOTE: This only works if tracks variable is
%                     called tracksFinal inside file.
%       minLength   : Minimum length of a track to be included in analysis.
%                     Optional. Default: 5 frames.
%       alpha       : Alpha-value for the statistical test to determine
%                     number of modes.
%                     Optional. Default: 0.01.
%       showPlot    : 0 to not plot anything.
%                     1 to plot the histogram and fitted exponentials.
%                     Optional. Default: 1.
%       maxNumMode  : Upper limit on the number of modes.
%                     Optional. Default: 10.            
%       binStrategy : Binning strategy for calculating the cumulative
%                     histogram. 1 for using "histogram" and 2 for using
%                     the data directly.
%                     Optional. Default: 2.
%       plotName    : The title of the plotted figure.
%                     Optional. Default: 'Figure'.
%
%OUTPUT modeParam   : Matrix with number of rows equal to number of modes
%                     and 4 columns:
%                     Column 1: diffusion coefficient of each mode
%                     Column 2: fraction of contribution of each mode
%                     Column 3: std of each diffusion coefficient
%                     Column 4: std of each fraction of contribution
%       expParam    : Output of fitHistWithExponentialsN
%       expControl  : Output of 
%
%REMARKS Code is written for 2D case only, but can be generalized to 3D.
%
%Khuloud Jaqaman, June 2012

%% Input

%get tracks if not input
%this only works if tracks are called tracksFinal inside the .mat file
if nargin < 1 || isempty(tracksFinal)
    [fName,dirName] = uigetfile('*.mat','Please choose .mat file of tracks.');
    load(fullfile(dirName,fName));
end

%check rest of input

if nargin < 2 || isempty(minLength)
    minLength = 5;
end

if nargin < 3 || isempty(alpha)
    alpha = 0.01;
end

if nargin < 4 || isempty(showPlot)
    showPlot = 1;
end

if nargin < 5 || isempty(maxNumMode)
    maxNumMode = 10;
end

if nargin < 6 || isempty(binStrategy)
    binStrategy = 2;
end

if nargin < 7 || isempty(plotName)
    plotName = 'Figure';
end

%% Mode decomposition

%keep only tracks of the right length
criteria.lifeTime.min = minLength;
indxKeep = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indxKeep);

%get number of tracks
numTracks = length(tracksFinal);

%get distribution of square displacements
%divide tracks into chunks of 5000 so as not to run out of memory
%also get localization precision
indxFirst = 1;
f2fdispSqAll = [];
xCoordStdAll = [];
yCoordStdAll = [];
while indxFirst < numTracks

    %get chunk of tracks for this iteration
    indxLast = min(indxFirst+5000-1,numTracks);
    tracksTmp = tracksFinal(indxFirst:indxLast);
    tracksMat = convStruct2MatIgnoreMS(tracksTmp,1);

    xCoord = tracksMat(:,1:8:end);
    yCoord = tracksMat(:,2:8:end);
    f2fdispSqTmp = diff(xCoord,1,2).^2 + diff(yCoord,1,2).^2;
    f2fdispSqTmp = f2fdispSqTmp(~isnan(f2fdispSqTmp));
    f2fdispSqAll = [f2fdispSqAll; f2fdispSqTmp]; %#ok<AGROW>

    xCoordStd = tracksMat(:,5:8:end);
    yCoordStd = tracksMat(:,6:8:end);
    xCoordStdAll = [xCoordStdAll; xCoordStd(~isnan(xCoordStd))]; %#ok<AGROW>
    yCoordStdAll = [yCoordStdAll; yCoordStd(~isnan(yCoordStd))]; %#ok<AGROW>
    
    indxFirst = indxLast + 1;

end
meanPosVar = mean([xCoordStdAll;yCoordStdAll].^2);

%fit exponentials to the distribution of square displacements
[~,binCenterP,expParam] = fitHistWithExponentialsN(f2fdispSqAll,alpha,showPlot,...
    maxNumMode,binStrategy,[],plotName,4*meanPosVar);

%calculate diffusion coefficient of each mode
%formula: mu (i.e. mean of exponential) = 4*diffCoef + 4*posStd^2
diffCoef = expParam(:,1)/4 - meanPosVar;

%calculate the diffusion coefficient std -- under-estimate
diffCoefStd = expParam(:,3)/4;

%calculate fraction of contribution of each mode
fracContr = expParam(:,2)/sum(expParam(:,2));

%calculate corresponding std -- under-estimate
fracContrStd = expParam(:,4)/sum(expParam(:,2));

%output
modeParam = [diffCoef fracContr diffCoefStd fracContrStd];

%test whether the distance between consecutive modes is significant given
%the diffusion coefficients and their standard deviations
if length(diffCoef) > 1
    diffCoefDelta = diff(diffCoef);
    deltaStd = sqrt( diffCoefStd(1:end-1).^2 + diffCoefStd(2:end).^2 );
    pvModeDiff = 1-normcdf(diffCoefDelta,0,deltaStd);
else
    pvModeDiff = NaN;
end

%control: repeat the above but for a synthetic data set generated from a
%mono-exponential distribution with the same average as the data-derived
%distribution of square displacements
numBin = length(binCenterP);
if binStrategy == 1
    binStrategyControl = 3;
else
    binStrategyControl = binStrategy;
end
f2fdispSqControl = exprnd(mean(f2fdispSqAll),length(f2fdispSqAll),1);
[~,~,expParamControl] = fitHistWithExponentialsN(f2fdispSqControl,alpha,showPlot,...
    maxNumMode,binStrategyControl,numBin,[plotName ' - Control'],4*meanPosVar);
diffCoefControl = expParamControl(:,1)/4 - meanPosVar;
diffCoefStdControl = expParamControl(:,3)/4;
fracContrControl = expParamControl(:,2)/sum(expParamControl(:,2));
fracContrStdControl = expParamControl(:,4)/sum(expParamControl(:,2));
modeParamControl = [diffCoefControl fracContrControl diffCoefStdControl fracContrStdControl];
if length(diffCoefControl) > 1
    diffCoefDeltaControl = diff(diffCoefControl);
    deltaStdControl = sqrt( diffCoefStdControl(1:end-1).^2 + diffCoefStdControl(2:end).^2 );
    pvModeDiffControl = 1-normcdf(diffCoefDeltaControl,0,deltaStdControl);
else
    pvModeDiffControl = NaN;
end

%% ~~~ the end ~~~



