function [Protrusion,Retraction] = getPersistenceTime(TS,deltaT,varargin)
%This function calculate the threshold to define a protrusion or retraction
%and uses it to estimate the time of protrusion and retraction for each
%event in TS
%Usage:
%       [Protrusion,Retraction,Up,Dw] = getPersistenceTime(TS,deltaT,varargin)
%
% Input:
%       TS     - vector - Time series
%       deltaT - scalar - sampling rate (in seconds, every deltaT seconds)
%       Optional:
%               per - number of standard deviations used to calculate the
%               threshold
%
% Output:
%       Protrusion.PersTime - vector - Protrusion time for all protrusive events
%       Protusion.BlockOut - cell   - Time point where Protrusion happened
%       RetrPersTime - vector - Retraction time for all retractive events
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
ip.addOptional('per',1,@isscalar);

ip.parse(TS,deltaT,varargin{:});
per  = ip.Results.per;
%**************************************************************************
%%

imf       = empiricalModeDecomp(TS)';
Mu        = mean(imf(1,:));
nPoint    = length(TS);
[~,noise] = testImf(imf);
sdtError  = std(noise); 

%Defining the lower and upper noise confidence bands
Protrusion.limit  = Mu + sdtError*per; 
Retraction.limit  = Mu - sdtError*per;
%*************************
 

TSprot  = NaN(size(TS));
TSretr  = NaN(size(TS));

TSprot(TS > Protrusion.limit) = TS(TS > Protrusion.limit);
TSretr(TS < Retraction.limit) = TS(TS < Retraction.limit);

ProtBlock = findBlock(setdiff(1:nPoint,find(isnan(TSprot))));
RetrBlock = findBlock(setdiff(1:nPoint,find(isnan(TSretr))));

Protrusion = getStuff(ProtBlock,TS,deltaT);
Retraction = getStuff(RetrBlock,-TS,deltaT);

end%End of main function

function cellData =  getStuff(block,TS,deltaT)

if ~isempty(ProtBlock)
    
    [cellData.PersTime,cellData.BlockOut,cellDa.MaxVeloc,cellDa.MeanVeloc] = findingProtRetrTime(block,TS,deltaT);
    
else
    
    cellDa.PersTime  = NaN;
    cellDa.BlockOut  = {[]};
    cellDa.MaxVeloc  = NaN;
    cellDa.MeanVeloc = NaN;
    
end

end