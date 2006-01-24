function [cutoffIndex, cutoffValue] = cutFirstHistMode(varargin);
%CUTFIRSTHISTMODE finds the end of the first mode in a histogram
%
% cutFirstHistMode is an implementation of the algorithm presented in "uni-
% modal thresholding" by P.L. Rosin, Pattern Recognition (2001); 34:2083.
% It assumes  that the first mode in the histogram (noise/background values
% for most applications) is strongest, and places the cutoff where distance
% between the line from the largest bin in the first mode to the first
% empty bin after the last nonempty bin and the histogram is largest.
%
% SYNOPSIS [cutoffIndex, cutoffValue] = cutFirstHistMode(counts, bins);
%          [cutoffIndex, cutoffValue] = cutFirstHistMode(data);
%          [...] = cutFirstHistMode(...,verbose)
%
% INPUT    counts, bins : counts in and center of histogram bins (output
%               of functions such as "hist" or "histogram". There has to be
%               more than one bin.
%          alternatively, you can pass the data directly, and the program
%               set up the histogram
%
%          verbose (optional) decides whether the function will open a
%               figure or not (default = 1).
%
%
% OUTPUT   cutoffIndex : Index into list of bins/data of the placement of
%               cutoff
%          cutoffValue : position of bin/data point where the histogram is
%               cut off
%
% REMARKS  If data is supplied, cutoffIndex/cutoffData will point to the
%               data point just at or above the center of the bin
%          The function works only on 1D data so far (multidimensional data
%          is handled as a vector)
%
%
% c: jonas, 8/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%======================
% TEST INPUT
%======================

% defaults
verbose = 1;

% goodDataIdx is used in case data contains nans
goodDataIdx = [];

switch nargin - isscalar(varargin{end})
    case 1 % data
        doHistogram = 1;
        data = varargin{1};
        data = data(:);

        if nargin == 2
            verbose = varargin{2};
        end
        
        % make sure that there are no nans or infs
        goodDataIdx = find(isfinite(data));
        if ~isempty(goodDataIdx)
            data = data(goodDataIdx);
        else
            error('all data is NaN or Inf!')
        end

    case 2 % conts,bins
        doHistogram = 0;
        counts = varargin{1};
        bins = varargin{2};
        counts = counts(:);
        bins = bins(:);
        nBins = length(bins);

        if nargin == 3
            verbose = varargin{3};
        end

    otherwise
        error('wrong number of input arguments or non-scalar ''verbose''')
end

%===========================


%==========================
% BUILD HISTOGRAM
%==========================
if doHistogram
    %     [counts,bins] = histogram(data);
    %     counts = counts(:);
    %     bins = bins(:);
    % instead of a histogram, make a cumulative histogram and flip it along the
    % horizontal. The number-count becomes counts, the actual value becomes the
    % bin position
    % we still need a histogram to find the maximum!

    %     [bins,sortIdx] = sort(data(:));
    %     [counts] = histogram(bins);
    %     [maxVal, maxIdx] = max(counts);
    %     bins = bins(maxIdx:end);
    %     nBins = length(bins);
    %     counts = [nBins:-1:1]';
    
    % the continuous histogram can potentially place the maximum elsewhere
    % (where it might belong, but that's not the point)
    % Do discrete histogram first, find max-bin, and then look for the
    % maximum in the continuous histogram only up to the end of the
    % maxBin+0.5
    [counts, bins] = histogram(data);
    % the last bin cannot become the maximum
    [dummy,maxBinIdx] = max(counts(1:end-1));
    maxBinValue = bins(maxBinIdx + 1);
    

    % continuous histogram
    [counts,bins] = histogram(data,'smooth');
    nBins = length(bins);
    
    % find maximum between bin 1 and the one with maxBinValue
    [maxVal, maxIdx] = max(counts(bins < maxBinValue));
    
else
    [maxVal, maxIdx] = max(counts);
end
%=========================


%=========================
% CUTOFF
%=========================

% for the line, we need the position of the maximum count and the position
% of the first empty bin following the last nonempty bin

% normalize counts and bins to go from 0 to 1.
countsN = counts / maxVal;
binsN = bins / max(abs(bins));

nBins = nBins - maxIdx + 1;
pointMax = [binsN(maxIdx), 1]; % maxVal = 1 now
pointEnd = [binsN(end) + median(diff(binsN)), countsN(end)];


% calculate perpendicular distance to the line
vector = pointEnd-pointMax;
[dummy,vector] = normList(vector);
distanceVector = perpVector(...
    repmat(pointMax,[nBins,1]),repmat(vector,[nBins,1]),...
    [binsN(maxIdx:end),countsN(maxIdx:end)]);
distance = normList(distanceVector);

% put distance = 0 wherever the contHistogram is above the line
distance(all(distanceVector>0,2)) = 0;

% find maximum
[maxDistance, maxDistanceIdx] = max(distance);

%========================


%========================
% ASSIGN OUTPUT
%========================

cutoffIndex = maxDistanceIdx + maxIdx - 1;
cutoffValue = bins(cutoffIndex);

if doHistogram

    % for plotting: stuff stolen from plotyy
    if verbose
        fig =figure;
        set(fig,'NextPlot','add')
        histogram(data)
        ax(1) = gca;
        ax(2) = axes('Units',get(ax(1),'Units'), ...
            'Position',get(ax(1),'Position'),'Parent',fig);
        plot(ax(2),bins,counts,'r',[cutoffValue,cutoffValue],[0,max(counts)],'-.r')
        set(ax(2),'YAxisLocation','right','Color','none','YColor','r')
        set(ax,'Box','off');
    end

else
    if verbose
        figure,bar(bins,counts),hold on,
        plot([cutoffValue,cutoffValue],[0,maxVal],'r')
    end

end

if ~isempty(goodDataIdx)
    cutoffIndex = goodDataIdx(cutoffIndex);
end