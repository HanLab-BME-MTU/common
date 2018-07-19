
function env = getEnvelope(In)
% Calculate the time series envelope 
% Used for EMD
% Marco Vilela, 2011

% Adding reflective boundary
In            = In(:);
maxLag        = 10;% Overshooting local neighborhood for a cubic spline
[acf,~,limit] = autocorr(In,maxLag);
bound         = max([max(find(acf>limit(1) | acf<limit(2)) ) 3]);
TS            = [flipud(In(2:bound));In;flipud(In(end - bound + 1:end - 1))];

npoint = length(TS);
[~,p]  = findpeaks(TS);

if numel(p) >= 2
    
    env =  spline(p,TS(p),1:npoint);
    env([1:bound-1 end-bound+2:end]) = [];
    
elseif ~isempty(p)
    
    npoint = length(In);
    [~,p]  = findpeaks(In);
    env    =  spline([0;p;npoint+1],[mean(In(1:2));In(p);mean(In(end-1:end))],1:npoint);
    
else
    
    env = 0;
    
end

end