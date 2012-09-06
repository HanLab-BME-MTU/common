function [ProtPersTime,RetrPersTime,ProtBlockOut,RetrBlockOut,Up,Dw] = getPersistenceTime(TS,deltaT,varargin)
%This function calculate the threshold to define a protrusion or retraction
%and uses it to estimate the time of protrusion and retraction for each
%event in TS
%Usage:
%       [ProtPersTime,RetrPersTime,ProtBlockOut,RetrBlockOut,Up,Dw] = getPersistenceTime(TS,deltaT,varargin)
%
% Input:
%       TS     - vector - Time series
%       deltaT - scalar - sampling rate (in seconds, every deltaT seconds)
%       Optional:
%               noiseStd - number of standard deviations used to calculate the
%                          local noise threshold (Default = 1)
%
% Output:
%       ProtPersTime - vector - Protrusion time for all protrusive events
%       RetrPersTime - vector - Retraction time for all retractive events
%       ProtBlockOut - cell   - Time point where Protrusion happened
%       RetrBlockOut - cell   - Time point where Retraction happened
%       Up           - scalar - Upper threshold
%       Dw           - scalar - Lower threshold
%
% See also: findingProtRetrTime, getEdgeMotionPersistence
%
%Marco Vilela, 2012

%% Parsing the input ******************************************************
ip = inputParser;
ip.addRequired('TS',@isvector);
ip.addRequired('deltaT',@isscalar);
ip.addOptional('noiseStd',1,@isscalar);
ip.parse(TS,deltaT,varargin{:});
noiseStd  = ip.Results.noiseStd;
%**************************************************************************
%%

imf       = emd(TS);
Mu        = mean(imf(1,:));
nPoint    = length(TS);
[~,noise] = testImf(imf);
sdtError  = std(noise); 

%Defining the lower and upper noise confidence bands
Up  = Mu + sdtError*noiseStd; 
Dw  = Mu - sdtError*noiseStd;
%*************************
 

TSprot = NaN(size(TS));
TSret  = NaN(size(TS));

TSprot(TS>Up) = TS(TS>Up);
TSret(TS<Dw)  = TS(TS<Dw);

ProtBlock = findBlock(setdiff(1:nPoint,find(isnan(TSprot))));
RetrBlock = findBlock(setdiff(1:nPoint,find(isnan(TSret))));
if ~isempty(ProtBlock)
    [ProtPersTime,ProtBlockOut] = findingProtRetrTime(ProtBlock,TS,deltaT);
else
    ProtPersTime = NaN;
    ProtBlockOut = {[]};
end
if ~isempty(RetrBlock)
    [RetrPersTime,RetrBlockOut] = findingProtRetrTime(RetrBlock,-TS,deltaT);
else
    RetrPersTime = NaN;
    RetrBlockOut = {[]};
end