function  [tau,tauAmp] = modifiedKendallCorr(x,varargin)
%Modified Kendall cross-correlation coefficient

ip = inputParser;
ip.addRequired('x',@(x) isvector(x));
ip.addOptional('y',x,@(x) isnumeric(x));
ip.addOptional('local',numel(x)-1,@isscalar);
ip.addOptional('alpha',0.05,@isscalar);
ip.addOptional('amp',false,@islogical);
ip.addOptional('maxLag',0,@isscalar);


ip.parse(x,varargin{:});
local    = ip.Results.local;
alpha    = ip.Results.alpha;
y        = ip.Results.y;
maxLag   = ip.Results.maxLag;

x = x(:);
y = y(:);

if ~all( size(x) == size(y) )
    error('x and y have different size')
else
    nObs = numel(x);
end

if (nObs - local) < maxLag
    disp('maximum lag is too large. Reducing it to the max allowed')
    maxLag = (nObs - local) - 1;
end


normalization       = local*(2*nObs - local - 1)/2;
contM               = getCorrMatrix(ones(nObs,1),local,[]);
[contX,ampM(:,:,1)] = getCorrMatrix(x,local,contM);
[contY,ampM(:,:,2)] = getCorrMatrix(y,local,contM);

%Matrix with [-1 0 1] for mismatch, no-match and match of slopes
match = contX.*contY;
tau   = sum(match(:))/normalization;

%Matrix with mean of [delta(i,j)] 
meanAmp  = mean(abs(ampM),3);
matchAmp = match.*meanAmp;
tauAmp   = sum(matchAmp(:))/sum(meanAmp(:));


contYr        = contY;
contXl        = contX;
ampM(:,:,3:4) = ampM(:,:,1:2);% 1 is for x and 2 for y
tauR          = NaN(1,maxLag-1);
tauL          = NaN(1,maxLag-1);
tauAmpR       = NaN(1,maxLag-1);
tauAmpL       = NaN(1,maxLag-1);


for iLag = 1:maxLag
    %Change this code - create a template that change with lag
    currNorm      = normalization - local*iLag;
    contYr        = circshift(contYr,-1);contYr(end,:) = 0;
    contXl        = circshift(contXl,-1);contXl(end,:) = 0;
    
    ampM(:,:,3)   = circshift(ampM(:,:,3),-1);ampM(end,:,3) = 0;
    ampM(:,:,4)   = circshift(ampM(:,:,4),-1);ampM(end,:,4) = 0;
    
    %Right - y series is shifted to the left
    matchR        = contX.*contYr;
    matchL        = contXl.*contY;
    
    
    meanAmpR      = mean(abs(ampM(:,:,[2 4])),3);
    matchAmpR     = matchR.*meanAmpR;
    tauAmpR(iLag) = sum(matchAmpR(:))/sum(meanAmpR(:));

    meanAmpL      = mean(abs(ampM(:,:,[1 3])),3);
    matchAmpL     = matchL.*meanAmpL;
    tauAmpL(iLag) = sum(matchAmpL(:))/sum(meanAmpL(:));
    
    %matchRamp     = matchR
    tauR(iLag)  = sum(matchR(:))/currNorm;
    tauL(iLag)  = sum(matchL(:))/currNorm;
    
end

tau    = [fliplr(tauL) tau tauR];
tauAmp = [fliplr(tauAmpL) tauAmp tauAmpR];
end%END OF MAIN FUNCTION

function [contTS,ampTSdiff] = getCorrMatrix(TS,neigh,contM)


hTS = hankel(TS);hTS(1,:) = [];
hTS(:,neigh+1:end)        = [];

if isempty(contM)   
    contTS    = hTS;
else
    ampTSdiff = (hTS - repmat(TS(1:end-1),1,neigh)).*contM;
    contTS    = sign(ampTSdiff);
end

end