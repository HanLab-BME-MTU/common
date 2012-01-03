function [avegCoh,CohCI,w]=coherenceBootstrap(S1,S2,varargin)
% This function calculates the power spectrum and coherence 
%
% Synopsis: 
%          []=coherenceBootstrap(S,S2,nWin,wType,noLap,Fs)
%
%Input:
%      S1    - signal
%      S2    - signal
%      nWin  - number of windows  
%      wType - window type. Options:
%              bartlett;
%              blackman;
%              blackmanharris;
%              bohmanwin;
%              chebwin;
%              flattopwin;
%              gausswin;
%              hamming; (Default)
%              hann;
%              kaiser;
%              nuttallwin;
%              parzenwin;
%              rectwin;
%              taylorwin;
%              triang;
%              tukeywin 
%      nOverLap - percentage of overlap between windows (Default 0.5 - 50%)
%
%Output

% Input check
ip=inputParser;
ip.addRequired('S1',@isnumeric);
ip.addRequired('S2',@(x) isnumeric(x) && isequal(size(x),size(S1)));
ip.addOptional('nWin',8,@isscalar);
ip.addOptional('wType','hamming',@ischar);
ip.addOptional('noLap',.5,@isscalar);
ip.addOptional('Fs',1,@isscalar);
ip.addParamValue('alpha',.05,@isscalar);
ip.addParamValue('nBoot',10000,@isscalar);
ip.parse(S1,S2,varargin{:})
nWin=ip.Results.nWin;
wType=ip.Results.wType;
noLap=ip.Results.noLap;
Fs=ip.Results.Fs;
alpha=ip.Results.alpha;
nBoot=ip.Results.nBoot;


[nPt,nVar] = size(S1);


q1         = alpha*100/2;
q2         = 100 - q1;
winLen     = floor( nPt/( (1 - noLap)*nWin + noLap ) );
exPt       = round(winLen*noLap);
windowF    = feval(wType,winLen);

%**************************************************************************
powerSpect = @(x) nanmean(abs(x).^2);
coherence  = @(x,y) abs( nanmean(y.*conj(x)) ).^2 ./ ( powerSpect(x).*powerSpect(y) );
%**************************************************************************

SS1 = nan(winLen,nWin);
SS2 = nan(winLen,nWin);

for i=1:nWin
    idx        = ((0:(nVar - 1))*nWin) + i;
    SS1(:,idx) = S1(1 + (winLen - exPt)*(i-1):(winLen - exPt)*i + exPt, :);
    SS2(:,idx) = S2(1 + (winLen - exPt)*(i-1):(winLen - exPt)*i + exPt, :);
end

%Tapering the signals
tSS1 = SS1.*repmat(windowF,1,nWin*nVar);
tSS2 = SS2.*repmat(windowF,1,nWin*nVar);

%FT the signals
nanIdx1 = find(sum(isnan(tSS1)));
nanIdx2 = find(sum(isnan(tSS2)));
nfft    = 2^nextpow2(nPt);
%Variables without NaN (Fast Fourier Transform)
dS1w(:,setdiff(1:nVar*nWin,nanIdx1)) = fft(tSS1(:, setdiff(1:nVar*nWin,nanIdx1) ),nfft)/nPt;
dS2w(:,setdiff(1:nVar*nWin,nanIdx2)) = fft(tSS2(:, setdiff(1:nVar*nWin,nanIdx2) ),nfft)/nPt;
%Variables with NaN (Extended Fourier Transform)
dS1w(:,nanIdx1) = edft(tSS1(:, nanIdx1 ),nfft)/nPt;
dS2w(:,nanIdx2) = edft(tSS2(:, nanIdx2 ),nfft)/nPt;

%One-sided FT
w            = Fs/2*linspace(0,1,nfft/2 +1);
oneSideSpec1 = 2*dS1w(1:nfft/2 +1,:);
oneSideSpec2 = 2*dS2w(1:nfft/2 +1,:);


opt = statset('UseParallel','always');
if nVar == 1
    avegCoh = feval(coherence,oneSideSpec1',oneSideSpec2');
    bootSp  = jackknife(coherence,oneSideSpec1',oneSideSpec2','Options',opt);
else
    
    bootSp  = tanh( bootstrp( nBoot, @(x,y) atanh(coherence(x,y)),oneSideSpec1',oneSideSpec2',...
                            'Options',opt) );
    avegCoh = nanmean(bootSp);
end
CohCI = prctile(bootSp,[q1 q2]);

