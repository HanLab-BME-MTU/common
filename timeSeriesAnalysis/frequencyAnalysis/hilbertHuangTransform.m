function [instAmp,instFreq,hht] = hilbertHuangTransform(TS,Fs,specRes)
%This function calculates the Hilbert spectrum using empirical
%mode decomposition as preliminary step
%
%USAGE:
%       [instAmp,instFreq,hht] = hilbertHuangTransform(TS,Fs,specRes)
%
%Input:
%       TS      - time series vector
%       Fs      - sampling frequency; real number in Hz
%       specRes - interger that specifies the number of frequency points in
%                 the final spectrum [frequency axis goes from 0 to Fs/2 with specRes number of points] 
%
%Output:
%       instAmp  - instantaneous Amplitute for each imf
%       instFreq - instantaneous frequency for each imf
%       hht      - structure: hht.timeVector [in seconds]
%                             hht.freqVector [in Hz]
%                             hht.spectrum  
%
%Ref:
%   book: Hilbert Transform and Application in Mechanical Vibration; Michael Feldman 
%
%Marco Vilela, 2013



%Decomposing signal into intrinsic modes
imf = emd(TS);

%Initialization
[nImf,nObs] = size(imf);
instAmp     = nan(nImf,nObs);
instFreq    = nan(nImf,nObs-1);
spectrum    = nan(specRes,nObs-1);

for i=1:nImf
    
    [instAmp(i,:),instFreq(i,:)] = getInstantaneousFrequency(imf(i,:));
    
end

instAmp(:,1) = [];%Make it same size as instFreq

%Guarantee frequencies within the physical range
instFreq = max(instFreq,0);
instFreq = min(instFreq,0.5);
% 0.5 is the Nyquist frequency for 1 unit of frequency. [No physical units for time]

%Constructing the time-frequency spectrum
dFreq                  = 0.5/(specRes - 1);
freqIdx                = round( instFreq./dFreq ) + 1;
finalFreqIdx           = sub2ind(size(spectrum),freqIdx,repmat(1:size(spectrum,2),size(freqIdx,1),1));
spectrum(finalFreqIdx) = instAmp;

%Adding unit to the instantaneous frequency [Hz]
instFreq       = Fs*instFreq;
hht.timeVector = (1/Fs)*(0:(nObs - 1));
hht.freqVector = Fs*(0:dFreq:0.5);
hht.spectrum   = spectrum;