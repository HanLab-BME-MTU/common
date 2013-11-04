function  [IA,IF] = getInstantaneousFrequency(TS)
%This function calculates the instantaneous amplitude and frequency of a
%monocromatic signal TS
%
%USAGE:
%       [IA,IF] = getInstantaneousFrequency(TS)
%
%Input:
%       TS - time series vector
%
%Output:
%       IA - instantaneous amplitude
%       IF - instantaneous frequency
%
% See also: hilbertHuangTransform
%
%Marco Vilela, 2013

%Analytical signal
analyticS = hilbert(TS);
IF        = abs(angle(analyticS(1:end-1).*conj(analyticS(2:end)))/(2*pi));
IA        = abs(analyticS);

end

