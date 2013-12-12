%readtiff() loads a tiff stack using libtiff
% This is ~1.5-2x faster than imread, useful for large stacks

% Francois Aguet, 05/21/2013

function s = readtiff(filepath, range, info)

w = warning('off', 'all'); % ignore unknown TIFF tags

if nargin<3 || isempty(info)
    info = imfinfo(filepath);
end

nx = info(1).Width;
ny = info(1).Height;

if nargin<2
    N = numel(info);
    range = 1:N;
else
    N = numel(range);
end

if info(1).BitDepth==16 && strcmpi(info(1).ColorType, 'grayscale')
    s = zeros(ny,nx,N,'uint16');
else
    s = zeros(ny,nx,N);
end

t = Tiff(filepath, 'r');
for i = 1:numel(range)
   t.setDirectory(range(i));
   s(:,:,i) = t.read();
end
t.close();

warning(w);
