function [sp,errStd,errV,varargout] = lsRandomFit(v1,v2,varargin)
%lsRandomFit: Spline (or linear) least-sqaure fitting between two random
%             variables with a large number of samples.
%
% Model: The goal of this function is to determine whether and how one random
%        variable depends on the other according to the following model:
%                         v2 = f(v1) + errV.
%        where v1 and errV approximately follow Gaussian distribution. The
%        outlier detection is based on a linear model (f being linear).
%        Therefore, if one wants to exclude outlier, f(v1) can not be too far
%        from being linear.
%
% SYNOPSIS:
%    [sp,errStd,errV] = lsRandomFit(v1,v2)
%    [sp,errStd,errV,outlierInd] = lsRandomFit(v1,v2,varargin)
%    [sp,errStd,errV,outlierInd,corrInd] = lsRandomFit(v1,v2,varargin)
%    [sp,errStd,errV,outlierInd,corrInd,exEllipse] = lsRandomFit(v1,v2,varargin)
%
% INPUT:
%    v1, v2: Vector of samples of two random variables.
%
%    Optional parameter/value pairs:
%       'outlierThreshold: Threshold for excluding outlier. It is a factor of 
%          the standard deviations of 'v1' and 'errV'.
%       'model': 'linear' or 'bspline'. When it is 'linear', a linear 
%                least-square fitting is calculated. The output is a line.
%       'figH'  : If a figure handle is given, the ellipse (for excluding
%                 outliers) and the fitting curve will be plotted. Outliers
%                 are also marked.
%
% OUTPUT:
%    sp : The sp-form of the fitting B-spline. If 'linear' is on, a line is 
%         given by sp = [a b] where v2 = a*v1+b+errV.
%    errStd: The standard deviation of the error vector 'errV'.
%
%    Optional output:
%       outlierInd: Indices of outliers.
%       exEllipse : Ellipse for excluding outliers. It is a structure that
%                   contains the following fields:
%                   eCenter, r1, r2, eAxis1, eAxis2, eAngle

if nargin < 2 | rem(nargin-2,2) ~= 0
   error('Wrong number of input arguments.');
end

outlierThreshold = Inf;
model            = 'bspline';
figH             = [];

if nargin > 2 
   for k = 1:2:nargin-2
      switch varargin{k}
         case 'outlierThreshold'
            outlierThreshold = varargin{k+1};
         case 'model'
            model = varargin{k+1};
         case 'figH'
            figH = varargin{k+1};
            if ~ishandle(figH)
               error('''figH'' is not a valid handle.');
            end
      end
   end
end

if ~strcmp(model,'linear') & ~strcmp(model,'bspline')
   error('The specified fitting model is not recognized.');
end

if isempty(outlierThreshold)
   outlierThreshold = Inf;
end

if isempty(v1) | isempty(v2)
    sp = [];
    errStd = [];
    errV   = [];
    outlierInd = [];
    corrInd    = [];
    
    r1 = NaN;
    r2 = NaN;
    eAxis1  = [];
    eAxis2  = [];
    eAngle  = NaN;
    eCenter = [];
    if nargout > 2
        varargout{1} = outlierInd;
    end

    if nargout > 3
        varargout{2} = corrInd;
    end

    if nargout > 4
        varargout{3}.center = eCenter;
        varargout{3}.r1      = r1;
        varargout{3}.r2      = r2;
        varargout{3}.axis1  = eAxis1;
        varargout{3}.axis2  = eAxis2;
        varargout{3}.angle  = eAngle;
    end
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

% %The minimum and maximum of 'v1' gives the range of fitting.
% minV1 = min(v1);
% maxV1 = max(v1);
% 
% %%%%%% First fit without excluding outlier
% numKnots = 4;
% errV = NaN*ones(size(v2));
% if strcmp(model,'bspline')
%    spOrder  = 4;
% 
%    %Create the knot sequence.
%    knots = augknt(linspace(minV1,maxV1,numKnots),spOrder);
% 
%    %Least sqaures spline fitting.
%    sp = spap2(knots,spOrder,v1, v2);
% 
%    errV = v2 - fnval(sp,v1);
% else
%    %Linear least-square fitting.
%    M = [v1.' ones(length(v1),1)];
%    sp = M.'*M\(M.'*v2.');
% 
%    errV = v2 - sp(1)*v1 - sp(2);
% end

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

corrInd = 1:length(v1);
for k = 1:numIterations
   %Calculate the covariance matrix.
   covM = cov(v1(corrInd),v2(corrInd));

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
   v1Mean = mean(v1(corrInd));
   v2Mean = mean(v2(corrInd));

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
      corrInd    = find(eDist<=2*outlierThreshold);
      outlierInd = find(eDist>2*outlierThreshold);
   elseif k == 2
      v1Iter2Std = r1;
      v2Iter2Std = r2;
      corrInd    = find(eDist<=outlierThreshold);
      outlierInd = find(eDist>outlierThreshold);
   else
      shrinkFactor = min(r1/v1Iter2Std,r2/v2Iter2Std);
      corrInd      = find(eDist<=outlierThreshold/shrinkFactor);
      outlierInd   = find(eDist>outlierThreshold/shrinkFactor);
      r1 = r1/shrinkFactor;
      r2 = r2/shrinkFactor;
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The minimum and maximum of 'v1' gives the range of fitting.
minV1 = min(v1(corrInd));
maxV1 = max(v1(corrInd));

%%%%%% Calculate the fitting curve. We use least-square B-spline.
numKnots = 4;
errV = NaN*ones(size(v2));
if strcmp(model,'bspline')
   spOrder  = 4;

   %Create the knot sequence.
   knots = augknt(linspace(minV1,maxV1,numKnots),spOrder);

   %Least sqaures spline fitting.
   sp = spap2(knots,spOrder,v1(corrInd), v2(corrInd));

   errV(corrInd) = v2(corrInd) - fnval(sp,v1(corrInd));
else
   %Linear least-square fitting.
   M = [v1(corrInd).' ones(length(v1(corrInd)),1)];
   sp = M.'*M\(M.'*v2(corrInd).');

   errV(corrInd) = v2(corrInd) - sp(1)*v1(corrInd) - sp(2);
end

errStd = std(errV(corrInd));
errV   = reshape(errV,size(v2));

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

   %Plot the fitting curve.
   plotRange = linspace(minV1,maxV1,5*numKnots);
   if strcmp(model,'bspline')
      plot(plotRange,fnval(sp,plotRange),'r','lineWidth',3);
   else
      plot(plotRange,sp(1)*plotRange+sp(2),'r','lineWidth',3);
   end
end

if nargout > 2
   varargout{1} = outlierInd;
end

if nargout > 3
   varargout{2} = corrInd;
end

if nargout > 4 
   varargout{3}.center = eCenter;
   varargout{3}.r1      = r1;
   varargout{3}.r2      = r2;
   varargout{3}.axis1  = eAxis1;
   varargout{3}.axis2  = eAxis2;
   varargout{3}.angle  = eAngle;
end
