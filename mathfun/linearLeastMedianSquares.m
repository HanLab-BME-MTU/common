function [x,Qxx,goodRows,sigmaB] = linearLeastMedianSquares(A,B,V,x0)
%LINEARLEASTMEDIANSQUARES uses a least median square estimate to do robust least squares fitting of data with outliers
%
% linearLeastMedianSquares fits linear problems in the form of A*x = B + E,
% where E is a hopefully gaussian error term. It first checks for outliers
% in the data (B), and discards those, then it fits the reduced data with a
% standard least squares algorithm (myLscov).
%
% see: see Danuser, 1992 or Rousseeuw & Leroy, 1987
%
% SYNOPSIS : [x,Qxx,goodRows,sigmaB] = linearLeastMedianSquares(A,B,V,x0)
%
% INPUT    : A,B : matrices to describe the least squares problem in the
%                  form A*x = B
%            V   : optional covariance matrix of the values in B (DO NOT INVERT)
%            x0  : optional initial guess (default: ones(size(x)))
%                  (fminsearch is sensitive to initial conditions, so it helps 
%                  to have a good first guess)
%
% OUTPUT   : x,Qxx : fitted unknowns and covariance
%            goodRows : rows in B that are not outliers
%            sigmaB   : estimate for the std of the error E without
%                       outliers
%            
% c: 3/04 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==============
% TEST INPUT
%==============

if nargin < 2 | isempty(A) | isempty(B)
    error('Not enough or empty input arguments in linearLeastMedianSquares');
end

% size A: necessary for defaults
sizA = size(A);

% V: default: eye
if nargin < 3 | isempty(V)
    V = eye(sizA(1));
    invV = V;
    diagInvV = diag(invV);
else
    % currently, no testing
    invV = inv(V);
    diagInvV = diag(invV);
end


% x0: default: ones
if nargin < 4 | isempty(x0)
    x0 = ones(sizA(2),1);
else
    % currently no testing
end

%===== END TEST INPUT ======


%========================
% LEAST MEDIAN SQUARES
%========================

% define magic numbers:
k=3; %value important for calculation of sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magicNumber2=1.4826^2; %see same publications

% generate function string
functionString = ['inline(''median((A*x-B).*diagInvV.*(A*x-B))'',''x'',''A'',''B'',''diagInvV'')'];
% generate rest of fminsearch call
fminCall = ['fminsearch(' functionString ',x0,options,A,B,diagInvV)'];

% set options
options = optimset('Display','off');

% minimize
xFmin = eval(fminCall);

% calculate statistics
res2 = (A*xFmin-B).^2;
medRes2 = median(res2);

%testvalue to calculate weights
testValue=res2/(magicNumber2*medRes2);

%goodRows: weight 1, badRows: weight 0
goodRows=find(testValue<=k^2);

% ssq=sum(res2);
sigmaB=sqrt(sum(res2(goodRows))/(length(goodRows)-4));

%====END LMS=========


%=======================
% LINEAR LEAST SQUARES
%=======================

% call myLscov
[x,Qxx,mse] = myLscov(A(goodRows,:),B(goodRows,:),V(goodRows,goodRows));