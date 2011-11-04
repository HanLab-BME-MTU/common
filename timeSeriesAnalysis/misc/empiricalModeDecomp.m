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
%       imf - intrinsic mode functions
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

[nvar,npoint] = size(X);

if nvar > npoint
    X = X';
    npoint = length(X);
end

%xTest   = [X(1) X X(end)];%That is because the matlab spline function sucks
xTest   = X;
imf     = [];
Snumber = [4 12];%The HHT and it's applications, Chapter 1, p,9
count1  = 1;

while ~isempty( findpeaks( xTest ) )
   x1      = xTest;
   ensImf  = zeros(Snumber(2),npoint);
   
   if sum( getEnvelope(-x1) ) ~= 0 && sum( getEnvelope(x1) ) ~= 0
       flag    = 0;
       %Creating a ensemble of IMF for each level
       count2 = 1;
       for i=1:Snumber(2)
           upperE = getEnvelope(x1);
           lowerE = -getEnvelope(-x1);
           x1     = x1 - nanmean([upperE;lowerE]);
           
           if imfTest( x1 )
               ensImf(i,:) = x1;
               flag        = 1;
               idx(count2) = i;
               count2      = count2 +1;
           end
           
       end
       
       if flag
           imf{end+1} = nanmean( ensImf( idx, : ) );
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



function out = imfTest(In)
%% Test if the input is a imf
N  = length(In);
%Number of zero-crossing
t1 = sum( In(1:N-1).*In(2:N) < 0);
%Numnber of extrema
t2  = length(findpeaks(In)) + length(findpeaks(-In));
out = 0;
if abs(t1 -t2) <= 1
    out = 1;
end


function env = getEnvelope(In)
%% calculate the time series envelope 
npoint = length(In);
[~,p]  = findpeaks(In);
%env    = spline([0 p npoint+1],[0 In(p) 0],1:npoint);
if ~isempty(p)
   % env    = spline([0 p npoint+1],[0 In(p) mean(In(end-1:end))],1:npoint);
    
   env =  spline([0 p npoint+1],[mean(In(1:2)) In(p) mean(In(end-1:end))],1:npoint);
else
    env = 0;
end
