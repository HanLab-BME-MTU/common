function  [tau,tauAmp] = modifiedKendallCorr(x,varargin)
%Modified Kendall cross-correlation coefficient

ip = inputParser;
ip.addRequired('x',@(x) isvector(x));
ip.addOptional('y',x,@(x) isnumeric(x));
ip.addOptional('neigh',numel(x)-1,@isscalar);
ip.addOptional('alpha',0.05,@isscalar);
ip.addOptional('amp',false,@islogical);
ip.addOptional('maxLag',5,@isscalar);


ip.parse(x,varargin{:});
neigh    = ip.Results.neigh;
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

if (nObs - neigh) < maxLag
    disp('maximum lag is too large. Reducing it to the max allowed')
    maxLag = (nObs - neigh) - 1;
end


normalization       = neigh*(2*nObs - neigh - 1)/2;
contM               = getCorrMatrix(ones(nObs,1),neigh,[]);
[contX,ampM(:,:,1)] = getCorrMatrix(x,neigh,contM);
[contY,ampM(:,:,2)] = getCorrMatrix(y,neigh,contM);

match = contX.*contY;
tau   = sum(match(:))/normalization;

matchAmp = match.*mean(abs(ampM),3);
tauAmp   = sum(matchAmp(:))/sum(abs(matchAmp(:)));


contYr = contY;
contXl = contX;
tauR   = NaN(1,maxLag-1);
tauL   = NaN(1,maxLag-1);

for iLag = 1:maxLag
    
    currNorm      = normalization - neigh*iLag;
    contYr        = circshift(contYr,-1);contYr(end,:) = 0;
    contXl        = circshift(contXl,-1);contXl(end,:) = 0;
    matchR        = contX.*contYr;
    matchL        = contXl.*contY;
    
    %matchRamp     = matchR
    tauR(iLag)  = sum(matchR(:))/currNorm;
    tauL(iLag)  = sum(matchL(:))/currNorm;
    
end

tau = [fliplr(tauL) tau tauR];

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