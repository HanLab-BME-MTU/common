function [tout,rout] = roseplot(varargin)
%ROSE   Angle histogram plot.
% Note: This function has been adapted from the Matlab 'rose' function and
% has been modified to include the following functionality: 
% 1) plot a normalized angle histogram. 
% 2) define radial axis limits.
%   
%   ROSE(THETA) plots the angle histogram for the angles in THETA.  
%   The angles in the vector THETA must be specified in radians.
%
%   ROSE(THETA,N) where N is a scalar, uses N equally spaced bins 
%   from 0 to 2*PI.  The default value for N is 20.
%
%   ROSE(THETA,X) where X is a vector, draws the histogram using the
%   bins specified in X.
%
%   ROSE(THETA,X,'scale','norm') normalizes the histogram to vary between
%   0 to 1.
%
%   ROSE(THETA,X,...,'rMax',rVal) sets the r-axis limit to rVal
%
%   ROSE(AX,...) plots into AX instead of GCA.
%
%   H = ROSE(...) returns a vector of line handles.
%
%   [T,R] = ROSE(...) returns the vectors T and R such that 
%   POLAR(T,R) is the histogram.  No plot is drawn.
%
%   See also HIST, POLAR, COMPASS.

%   Clay M. Thompson 7-9-91
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.14.4.4 $  $Date: 2005/04/28 19:56:53 $

[cax,args,nargs] = axescheck(varargin{:});    
error(nargchk(1,6,nargs,'struct'));

theta = args{1};
if nargs > 1, 
  x = args{2}; 
end

if ischar(theta)
  error(id('NonNumericInput'),'Input arguments must be numeric.');
end
theta = rem(rem(theta,2*pi)+2*pi,2*pi); % Make sure 0 <= theta <= 2*pi
if nargs==1,
  x = (0:19)*pi/10+pi/20;

elseif nargs>=2,
  if ischar(x)
    error(id('NonNumericInput'),'Input arguments must be numeric.');
  end
  if length(x)==1,
    x = (0:x-1)*2*pi/x + pi/x;
  else
    x = sort(rem(x(:)',2*pi));
  end
end
if nargs>=4 && strcmpi(args{3},'scale') && strcmpi(args{4},'norm')
    doNorm=1;
else
    doNorm=0;
end

if nargs>=6 && strcmpi(args{5},'rMax') && isnumeric(args{6}) && ~isempty(args{6})
    setRmax=1;
    rMax=args{6};
else
    setRmax=0;
end
    
if ischar(x) || ischar(theta)
  error(id('NonNumericInput'),'Input arguments must be numeric.');
end

% Determine bin edges and get histogram
edges = sort(rem([(x(2:end)+x(1:end-1))/2 (x(end)+x(1)+2*pi)/2],2*pi));
edges = [edges edges(1)+2*pi];
nn = histc(rem(theta+2*pi-edges(1),2*pi),edges-edges(1));
nn(end-1) = nn(end-1)+nn(end);
nn(end) = [];

% Form radius values for histogram triangle
if min(size(nn))==1, % Vector
  nn = nn(:); 
end
[m,n] = size(nn);
mm = 4*m;
r = zeros(mm,n);
r(2:4:mm,:) = nn;
r(3:4:mm,:) = nn;

if doNorm
    r=r/sum(nn);
end

% Form theta values for histogram triangle from triangle centers (xx)
zz = edges;

t = zeros(mm,1);
t(2:4:mm) = zz(1:m);
t(3:4:mm) = zz(2:m+1);

if nargout<2
  % Generate the correct axis limits:
  if setRmax
      h_fake = polar(0,rMax);
      hold on;
  end
  
  if ~isempty(cax)
    h = polar(cax,t,r);
  else
    h = polar(t,r);
  end
  
  if setRmax
    set(h_fake, 'Visible', 'Off');
  end

  
  % If someone knows how to fix this please do it:
  % Register handles with m-code generator
  % if ~isempty(h)
  %    mcoderegister('Handles',h,'Target',h(1),'Name','rose');
  % end
  
  if nargout==1, tout = h; end
  return
end

if min(size(nn))==1,
  tout = t'; rout = r';
else
  tout = t; rout = r;
end

function str=id(str)
str = ['MATLAB:rose:' str];