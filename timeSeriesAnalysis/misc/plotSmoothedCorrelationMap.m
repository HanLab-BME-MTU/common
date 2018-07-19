function  plotSmoothedCorrelationMap(xCorr,Fs)
% Plots smoothed version of the input xCorr
%
%Usage:
%       plotSmoothedCorrelationMap(xCorr,Fs)
%
%Input:
%       xCorr - correlation map - (2*#Lags+1,#Variables)
%       Fs    - sampling rate 
%
%Marco Vilela, 2013

nLag   = size(xCorr,1);
maxLag = (nLag-1)/2;
sCC    = smoothActivityMap(xCorr);

figure
imagesc((1:size(sCC,2)/5),Fs.*linspace(-maxLag,maxLag,size(sCC,1)),sCC)


ylabel('Lag','FontSize',14)
xlabel('Windows','FontSize',14)
