function [aAcf,pAcf,CC]=windowsTimeSeriesTest(movieData,inCell,baTch)
%This function runs calculates the autocorrelation and cross-correlation
%between the protrusion and activity maps
%Input :
%       movieData - movieData object
%       inCell    - layer of windows to be used as activity Ex: inCell = 1
%       means cell edge.Default = 1;
%       baTch     - Suppresses the plot if set to 1. Default = 0;
%Output:
%       aAcf - Activity map autocorrelation
%       pAcf - Protrusion map autocorrelation
%       CC   - CrossCorrelation between protrusion and activity
%Marco Vilela, 2011

if ~isa(movieData,'MovieData')
    error('Input is not a moviedata object');
end
if nargin <2
    inCell = 1;
end
if nargin < 3
    baTch = 0;
end

numProc=length(movieData.processes_);
idxAct=[];
idxProt=[];

for i=1:numProc
    if isa(movieData.processes_{i},'ProtrusionSamplingProcess');
        idxProt=i;
    elseif isa(movieData.processes_{i},'WindowSamplingProcess');
        idxAct=i;
    end
    if ~isempty(idxProt) && ~isempty(idxAct)
        break;
    end
end

%Load cell protrusion and protein activity
Prot       = movieData.processes_{idxProt}.loadChannelOutput;
protrusion = Prot.avgNormal;
Act        = movieData.processes_{idxAct}.loadChannelOutput(1);
activity   = squeeze(Act.avg(:,inCell,:));%Gets the activity at the inCell layer

[nWind,nPoint] = size(protrusion);

%% Outliers detection
protrusion(detectOutliers(protrusion,5)) = NaN;
activity(detectOutliers(activity,5))     = NaN;

%% Percentage of NaN
[Prot,Prange] = removeMeanTrendNaN(protrusion');
[Act,Arange]  = removeMeanTrendNaN(activity');

%% Autocorrelation
%At least 50 points are needed to calculate the ACF
%Number of lags <= N/4;
%Ref: Time Series Analysis, Forecast and Control. Jenkins, G. Box,G
aAcf     = nan(round(nPoint/4),nWind);
pAcf     = nan(round(nPoint/4),nWind);
CC       = nan(round(nPoint/4),nWind);
abounds  = nan(2,nWind);
pbounds  = nan(2,nWind);
ccBounds = nan(2,nWind);
minP     = 50;
for i=1:nWind
    if length(Act{i}) >= minP
        [aAcf(1:round(length(Act{i})/4)+1,i),~,abounds(:,i)] = autocorr(Act{i},round(length(Act{i})/4));
    end
    if length(Prot{i}) >= minP
        [pAcf(1:round(length(Prot{i})/4)+1,i),~,pbounds(:,i)] = autocorr(Prot{i},round(length(Prot{i})/4));
    end
    
    [~,rangeP,rangeA] = intersect(Prange{i},Arange{i});
    ccL               = length(rangeP);
    if ccL >= minP
        [CC(1:2*round( ccL/4 )+1,i),~,ccBounds(:,i)] = crosscorr(Prot{i}(rangeP),Act{i}(rangeA),round( ccL/4 ) );
    end
end

if ~baTch
%% Plot Activity autocorrelation    
    mesh(aAcf,'FaceColor','interp');hold on
    [acL,~] = size(aAcf);
    upline  = repmat(abounds(1,:),acL,1);
    mesh(upline,'FaceColor',[63/255 63/255 63/255])
    
    dline  = repmat(abounds(2,:),acL,1);
    mesh(dline,'FaceColor',[63/255 63/255 63/255])
    
    ylabel('Lag')
    xlabel('Window')
    zlabel('Activity Autocorrelation')
    
    clear upline;clear dline
%% Plot Protrusion autocorrelation   
    figure
    mesh(pAcf,'FaceColor','interp');hold on
    [pcL,~] = size(pAcf);
    
    upline  = repmat(pbounds(1,:),pcL,1);
    mesh(upline,'FaceColor',[63/255 63/255 63/255])
    
    dline  = repmat(pbounds(2,:),pcL,1);
    mesh(dline,'FaceColor',[63/255 63/255 63/255])
    
    ylabel('Lag')
    xlabel('Window')
    zlabel('Protrusion Autocorrelation')
    
    clear upline;clear dline
  %% Plot CrossCorrelation 
    figure
    mesh(CC,'FaceColor','interp');hold on
    [cL,~] = size(CC);
    
    upline  = repmat(ccBounds(1,:),cL,1);
    mesh(upline,'FaceColor',[63/255 63/255 63/255])
    
    dline  = repmat(ccBounds(2,:),cL,1);
    mesh(dline,'FaceColor',[63/255 63/255 63/255])
    
    ylabel('Lag')
    xlabel('Window')
    zlabel('Crosscorrelation(Protrusion,Activity)')
end