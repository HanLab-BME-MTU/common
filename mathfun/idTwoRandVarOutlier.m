function [outlier,inlier,exEllipse] = idTwoRandVarOutlier(v1,v2,varargin)
%idTwoRandVarOutlier: Identify outlier samples of two random variables.
%
% Note: This function uses the eigenvalues and eigenvectors of the
% covariance matrix of two random variables to draw an ellipse of threshod 
% radii. Population inside the ellipse is identified as the correlated
% sampling cloud of the two random variables.
%
% SYNOPSIS:
%    [outlier,inlier,exEllipse] = lsRandomFit(v1,v2,varargin)
%
% INPUT:
%    v1, v2: Vector of samples of two random variables.
%
%    Optional parameter/value pairs:
%       'threshold' : Threshold for excluding outlier. It is a multiplying 
%                     factor of the standard deviations of 'v1' and 'v2'.
%                     Default: 4.
%       'figH'       : If a figure handle is given, the ellipse (for excluding
%                      outliers) will be plotted. Outliers are also marked.
%
% OUTPUT:
%       outlier   : Indices of outliers.
%       inlier    : Indices of correlated samples inside the ellipse.
%       exEllipse : Ellipse for excluding outliers. It is a structure that
%                   contains the following fields:
%                   eCenter : The center of the ellipse.
%                   r1, r2  : The two radii of the ellipse.
%                   eAxis1,
%                   eAxis2  : The two axis of the ellipse.
%                   eAngle  : The angle between the long axis and x-axis.

if nargin < 2 | rem(nargin-2,2) ~= 0
   error('Wrong number of input arguments.');
end

outlierThreshold = 4;
figH             = [];

if nargin > 2 
   for k = 1:2:nargin-2
      switch varargin{k}
         case 'threshold'
            outlierThreshold = varargin{k+1};
         case 'figH'
            figH = varargin{k+1};
            if ~ishandle(figH)
               error('''figH'' is not a valid handle.');
            end
      end
   end
end

if isempty(v1) | isempty(v2)
    outlier = [];
    inlier  = [];

    exEllipse = [];
    return;
end

[m,n] = size(v1);
if n == 1
   v1 = v1.';
end
[m,n] = size(v2);
if n == 1
   v2 = v2.';
end


%%%%%%% Exclude outlier
%Assumming the 'v1' and 'v2' are linear transformation of two independent
% Gaussian distribution, then the eigenvalues of the covariance matrix of
% 'v1' and 'v2' gives the variances of the two independent Gaussian
% distributions. Two independent random variables of Gaussian distribution
% form an elliptic shape on the 2D plane when they are sampled with a large
% population. The ratio between the long and short radius is given by the
% ratio between the two standard deviations of the two Guassian distributions. 
% The bigger the ellipse, the lower the probability of being
% outside the ellipse. We then exclude those (v1,v2) samples that are outside 
% an ellipse determined by a threshold factor of the standard deviations of 
% the two independent random variables.
%

numIterations = 2;

inlier = 1:length(v1);
for k = 1:numIterations
   %Calculate the covariance matrix.
   covM = cov(v1(inlier),v2(inlier));

   % 'S' is the linear transformation.
   [S,D] = eig(covM);

   %The variances of the two independent random variables.
   if D(1,1) < D(2,2)
      %Make sure the first eigenvalue is alway bigger.
      tmp = D(1,1);
      D(1,1) = D(2,2);
      D(2,2) = tmp;

      tmp = S(:,1);
      S(:,1) = S(:,2);
      S(:,2) = tmp;
   end

   var1 = D(1,1);
   var2 = D(2,2);

   %The two axes of the ellipse.
   eAxis1 = S(:,1);
   eAxis2 = S(:,2);

   %The two radii of the threshold ellipse for excluding outliers.
   r1 = sqrt(var1);
   r2 = sqrt(var2);

   %The mean of the two input vectors.
   v1Mean = mean(v1(inlier));
   v2Mean = mean(v2(inlier));

   %Center of the ellipse.
   eCenter = [v1Mean v2Mean].';

   %Angle of the long axis.
   eAngle = acos(sign(eAxis1(2))*eAxis1(1)/norm(eAxis1));

   %Calculate elliptic distance of each (v1,v2) data point to the center of 
   % the ellipse. First, map all the (v1,v2) to the two independent random 
   % variables.
   randV = S.'*([v1-v1Mean; v2-v2Mean]);
   eDist = sqrt(randV(1,:).^2/var1 + randV(2,:).^2/var2);

   if k == 1
      inlier  = find(eDist<=2*outlierThreshold);
      outlier = find(eDist>2*outlierThreshold);
   elseif k == 2
      v1Iter2Std = r1;
      v2Iter2Std = r2;
      inlier  = find(eDist<=outlierThreshold);
      outlier = find(eDist>outlierThreshold);
   else
      shrinkFactor = min(r1/v1Iter2Std,r2/v2Iter2Std);
      inlier  = find(eDist<=outlierThreshold/shrinkFactor);
      outlier = find(eDist>outlierThreshold/shrinkFactor);
      r1 = r1/shrinkFactor;
      r2 = r2/shrinkFactor;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ishandle(figH)
   figure(figH); hold off;

   plot(v1(corrInd),v2(corrInd),'.'); hold on;
   plot(v1(outlierInd),v2(outlierInd),'r.');

   plotellipse(eCenter,outlierThreshold*r1,outlierThreshold*r2,eAngle,'m');

   %Plot the two axes
   plot(eCenter(1)+outlierThreshold*r1*eAxis1(1)*[1 -1], ...
      eCenter(2)+outlierThreshold*r1*eAxis1(2)*[1 -1],'m');
   plot(eCenter(1)+outlierThreshold*r2*eAxis2(1)*[1 -1], ...
      eCenter(2)+outlierThreshold*r2*eAxis2(2)*[1 -1],'m');
end

exEllipse.center = eCenter;
exEllipse.r1     = r1;
exEllipse.r2     = r2;
exEllipse.axis1  = eAxis1;
exEllipse.axis2  = eAxis2;
exEllipse.angle  = eAngle;
