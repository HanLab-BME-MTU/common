function [imf,flag] = empiricalModeDecomp(X)
% This function calculates a set of intrinsic mode function from the input
% signal X
%
%Synopsis:
%         imf = empiricalModeDecomp(X)  
%
% Input:
%       X - 1-D column vector - time series
%Output
%       imf  - intrinsic mode functions
%       flag - 1 when is Ok; 0 when the spline is fucked up and no imf is
%              found
%
%References:
%
%N. E. Huang, Z. Shen, and S. R. Long et al. The empirical mode
%decomposition and the Hilbert spectrum for nonlinear and non-stationary
%time series analysis.
%Proceedings of the Royal Society of London,
%A(454):903?995, 1998.
%
%Marco Vilela, 2011
X       = X(:);
nObs    = numel(X);
xTest   = X;
imf     = [];
Snumber = [4 12];%The HHT and it's applications, Chapter 1, p,9
count1  = 1;
flag    = [];
while ~isempty( findpeaks( xTest ) )
   x1      = xTest;
   ensImf  = zeros(nObs,Snumber(2));
   
   if sum( getEnvelope(-x1) ) ~= 0 && sum( getEnvelope(x1) ) ~= 0
       flag    = 0;
       %Creating a ensemble of IMF for each level
       count2 = 1;
       for i=1:Snumber(2)
           upperE = getEnvelope(x1)';
           lowerE = -getEnvelope(-x1)';
           x1     = x1 - nanmean([upperE lowerE],2);
           
           if imfTest( x1 )
               
               ensImf(:,count2) = x1;
               flag             = 1;
               count2           = count2 +1;
               
           end
           
       end
       
       if flag
           imf{end+1} = nanmean( ensImf, 2 );
       else
           break;
       end
       
   else
       imf{end+1} = x1;
   end
   
   xTest  = xTest - imf{count1};
   count1 = count1 + 1;
   clear ensImf;
end

imf{end+1} = xTest;

imf = cell2mat(imf);

end%End of main function


%% Test if the input is a imf
function out = imfTest(In)

N  = length(In);

%Number of zero-crossing
t1 = sum( In(1:N-1).*In(2:N) < 0);

%Numnber of extrema
t2  = length(findpeaks(In,'minpeakdistance',2)) + length(findpeaks(-In,'minpeakdistance',2));
out = 0;

if abs(t1 -t2) <= 1
    
    out = 1;
    
end

end%End of imfTest

%% Calculate the time series envelope 
function env = getEnvelope(In)

%Adding reflective boundary
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

end%End of getEnvelope