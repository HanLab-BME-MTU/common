function mergeFrames(path1, path2, mergepath)

if nargin<1
    path1 = uigetdir('.tif', 'Select directory frames of 1st movie');
end
if nargin<2
    path2 = uigetdir('.tif', 'Select directory frames of 2nd movie');
end
if nargin<3
    mergepath = uigetdir('.tif', 'Select directory for output');
end

disp('adw')

tifFiles1 = dir([path1 '*.tif*']);
tifFiles2 = dir([path2 '*.tif*']);

nFrames = length(tifFiles1);

for k = 1:nFrames
    img1 = double(imread([path1 tifFiles1(k).name]));
    img2 = double(imread([path2 tifFiles2(k).name]));
    [nx ny] = size(img1);
    merge = zeros(nx, 2*ny);
    merge(:,1:ny) = img1;
    merge(:,ny+1:2*ny) = img2;
    figure; imagesc(merge)
    imwrite(uint16(merge), [mergepath 'mergeframe_' num2str(k, ['%.' num2str(length(num2str(nFrames))) 'd']) '.tif'], 'tif', 'compression', 'lzw');
end