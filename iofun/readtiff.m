%readtiff() loads a tiff stack using libtiff
% This is ~1.5-2x faster than imread, useful for large stacks

% Francois Aguet, 05/21/2013

function s = readtiff(filepath, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('filepath');
ip.addOptional('range', []);
ip.addOptional('info', [], @isstruct);
ip.addParamValue('ShowWaitbar', false, @islogical);
ip.parse(filepath, varargin{:});
info = ip.Results.info;
range = ip.Results.range;
if isempty(info)
    info = imfinfo(filepath);
end
if isempty(range)
    N = numel(info);
    range = 1:N;
else
    N = numel(range);
end

w = warning('off', 'all'); % ignore unknown TIFF tags

nx = info(1).Width;
ny = info(1).Height;

if info(1).BitDepth==16 && strcmpi(info(1).ColorType, 'grayscale')
    s = zeros(ny,nx,N,'uint16');
else
    s = zeros(ny,nx,N);
end

t = Tiff(filepath, 'r');
if ~ip.Results.ShowWaitbar
    for i = 1:numel(range)
        t.setDirectory(range(i));
        s(:,:,i) = t.read();
    end
else
    [~,fname] = fileparts(filepath);
    h = waitbar(0, ['Loading ' fname]);
    for i = 1:numel(range)
        t.setDirectory(range(i));
        s(:,:,i) = t.read();
        waitbar(i/numel(range))
    end
    close(h);
end
t.close();

warning(w);
