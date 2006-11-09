function [N,X,sp] = histogram(varargin)
% HISTOGRAM generates a histogram using the "optimal" number of bins
%
% If called with no output argument, histogram plots into the current axes
%
% SYNOPSIS [N,X,sp] = histogram(data,factor)
%          [...] = histogram(data,'smooth')
%          [...] = histogram(axesHandle,...)
%
% INPUT    data: vector of input data
%          factor: (opt) factor by which the bin-widths are multiplied
%                   if "smooth", a smooth histogram will be formed.
%                   (requires the spline toolbox)
%          axesHandle: (opt) if given, histogram will be plotted into these
%                       axes, even if output arguments are requested
%
% OUTPUT   N   : number of points per bin (value of spline)
%          X   : center position of bins (sorted input data)
%          sp  : definition of the smooth spline
%
% REMARKS: The smooth histogram is formed by calculating the cumulative
%           histogram, fitting it with a smoothening spline and then taking
%           the analytical derivative. If the number of data points is
%           markedly above 1000, the spline is fitting the curve too
%           locally, so that the derivative can have huge peaks. Therefore,
%           only 1000-1999 points are used for estimation.
%           Note that the integral of the spline is almost exactly the
%           total number of data points. For a standard histogram, the sum
%           of the hights of the bins (but not their integral) equals the
%           total number of data points. Therefore, the counts might seem
%           off.
%
%           WARNING: If there are multiples of the minimum value, the
%           smooth histogram might get very steep at the beginning and
%           produce an unwanted peak. In such a case, remove the
%           multiple small values first (for example, using isApproxEqual)

%
% c: 2/05 jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
if nargin < 1
    error('not enough input arguments for histogram')
end

% check for axes handle
if length(varargin{1}) == 1 && ishandle(varargin{1});
    axesHandle = varargin{1};
    varargin(1) = [];
else
    % ensure compatibility to when axesHandle was given as last input
    if nargin == 3 && ishandle(varargin{end})
        axesHandle = varargin{end};
        varargin(end) = [];
    else
        axesHandle = 0;
    end
end

% assign data
numArgIn = length(varargin);
data = varargin{1};
data = data(:);

% check for non-finite data points
data(~isfinite(data)) = [];

% check for "factor"
if numArgIn < 2 || isempty(varargin{2})
    factor = 1;
else
    factor = varargin{2};
end
if ischar(factor)
    if strmatch(factor,'smooth')
        factor = -1;
    else
        error('only string input permitted is ''smooth''')
    end
end

% doPlot is set to 1 for now. We change it to 0 below if necessary.
doPlot = 1;

nData = length(data);
% check whether we do a standard or a smooth histogram
if factor ~= -1
    if nData < 20
        warning('HISTOGRAM:notEnoughDataPoints','Less than 20 data points!')
        nBins = ceil(nData/4);
    else


        % create bins with the optimal bin width
        % W = 2*(IQD)*N^(-1/3)
        interQuartileDist = diff(prctile(data,[25,75]));
        binLength = 2*interQuartileDist*length(data)^(-1/3)*factor;

        % number of bins: divide data range by binLength
        nBins = floor((max(data)-min(data))/binLength);

        if ~isfinite(nBins)
            nBins = length(unique(data));
        end

    end

    % histogram
    if nargout > 0
        [N,X] = hist(data,nBins);
        % adjust the height of the histogram
        Z = trapz(X,N);
        N = N * nData/Z;
    else
        hist(data,nBins);
    end

else
    % make cdf, smooth with spline, then take the derivative of the spline

    % cdf
    xData = sort(data);
    yData = 1:nData;

    % when using too many data points, the spline fits very locally, and
    % the derivatives can still be huge. Good results can be obtained with
    % 500-1000 points. Use 1000 for now
    step = max(floor(nData/1000),1);
    xData2 = xData(1:step:end);
    yData2 = yData(1:step:end);

    % spline. Use strong smoothing
    cdfSpline = csaps(xData2,yData2,1./(1+mean(diff(xData2))^3/0.0006));

    % pdf is the derivative of the cdf
    pdfSpline = fnder(cdfSpline);

    % histogram
    if nargout > 0
        xDataU = unique(xData);
        N = fnval(pdfSpline,xDataU);
        X = xDataU;
        % adjust the height of the histogram
        Z = trapz(X,N);
        N = N * nData/Z;
        sp = pdfSpline;
        % set doPlot. If there is an axesHandle, we will plot
        doPlot = axesHandle;
    end
    % check if we have to plot. If we assigned an output, there will only
    % be plotting if there is an axesHandle.
    if doPlot
        if axesHandle
            plot(axesHandle,xData,fnval(pdfSpline,xData));
        else
            plot(xData,fnval(pdfSpline,xData));
        end
    end
end