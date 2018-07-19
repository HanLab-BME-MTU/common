function hh = herrorbar(varargin)
%HERRORBAR Horizontal error bar plot.
% This function draws horizontal error bars by calling ERRORBAR and rotating the result.
% See ERRORBAR for usage instructions/documentation.
%
% Additional option: 'BarWidth': half-size of the errorbar 'T'
% 
% Example: 
% herrorbar(1:3, [3 2 4], [0.5 1 1.5], 'r', 'LineStyle', 'none', 'BarWidth', 0.3);

% Francois Aguet, 01/21/2012.

idx = find(strcmpi(varargin, 'BarWidth'));
if ~isempty(idx)
    de = varargin{idx+1};
    varargin([idx idx+1]) = [];
    hh = errorbar(varargin{:});
    setErrorbarStyle(hh, de);
else
    hh = errorbar(varargin{:});
end

hc = get(hh, 'Children');

N = numel(get(hc(1), 'XData'));

if nargin<3
    x = 1:N;
    y = varargin{1};
else
    x = varargin{1};
    y = varargin{2};
end

xv = reshape(repmat(x, [9 1]), [1 9*N]);
yv = reshape(repmat(y, [9 1]), [1 9*N]);
xdata = get(hc(2), 'YData') - yv + xv;
ydata = get(hc(2), 'XData') - xv + yv;

set(hc(2), 'XData', xdata, 'YData', ydata);