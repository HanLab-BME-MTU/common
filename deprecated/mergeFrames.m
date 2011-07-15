function mergeFrames(path1, path2, mergepath)

fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

if nargin<1
    path1 = uigetdir('Select directory frames of 1st movie');
end
if nargin<2
    path2 = uigetdir('Select directory frames of 2nd movie');
end
if nargin<3
    mergepath = uigetdir('Select directory for output');
end

tifFiles1 = imDir([path1 filesep]);
tifFiles2 = imDir([path2 filesep]);

nFrames = length(tifFiles1);

for k = 1:nFrames
    img1 = double(imread([path1 filesep tifFiles1(k).name]));
    img2 = double(imread([path2 filesep tifFiles2(k).name]));
    [nx ny] = size(img1);
    merge = zeros(nx, 2*ny);
    merge(:,1:ny) = img1;
    merge(:,ny+1:2*ny) = img2;
    imwrite(uint16(merge), [mergepath filesep 'mergeframe_' num2str(k, ['%.' num2str(length(num2str(nFrames))) 'd']) '.tif'], 'tif', 'compression', 'lzw');
end