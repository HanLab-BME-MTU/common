function mask = medianSeg(im, closureRadius)

img_org = im;
h = fspecial('gaussian',[7 7], 2);
img_org = imfilter(img_org,h,'replicate');

imf = im(:);
Q(2) = median(imf(:));
imf2(find(imf > Q(2))) = 0;
img_org(find(im>Q(2))) = Q(2);

thismask = nrm(img_org,16);
level = graythresh(thismask);

if level > 0
    mask = im2bw(thismask,level);
    L = bwlabel(mask);
    s  = regionprops(L, 'Area');
    Allarea = [s.Area];
    [tmp1, tmp2] = max(Allarea);
    mask = (L == tmp2);
    
    closureBrush = strel('disk',closureRadius);
    mask = imclose(cast(mask,'double'),closureBrush);
    mask = cast(mask,'logical');

    mask = double(imfill(double(mask),'holes'));
    
    if length(unique(mask)) == 1
        mask = thismask > median(double(thismask(:)));
        L = bwlabel(mask);
        s  = regionprops(L, 'Area');
        Allarea = [s.Area];
        [tmp1, tmp2] = max(Allarea);
        mask = (L == tmp2);

        closureBrush = strel('disk',closureRadius);
        mask = imclose(cast(mask,'double'),closureBrush);
        mask = cast(mask,'logical');

        mask = double(imfill(double(mask),'holes'));

    end
else
    mask = thismask > median(thismask(:));
    L = bwlabel(mask);
    s  = regionprops(L, 'Area');
    Allarea = [s.Area];
    [tmp1, tmp2] = max(Allarea);
    mask = (L == tmp2);

    closureBrush = strel('disk',closureRadius);
    mask = imclose(cast(mask,'double'),closureBrush);
    mask = cast(mask,'logical');

    mask = double(imfill(double(mask),'holes'));
end

L2 = bwlabel(mask==0);
s2  = regionprops(L2, 'Area','PixelIdxList');

numberOfBlobs = size(s2, 1);
for k=1:numberOfBlobs
    PixelIdxList = s2(k).PixelIdxList;
    meanGL(k) = mean(thismask(PixelIdxList));
end
[tmp1, tmp2] = min(meanGL);
maskbg = (L2 == tmp2);
mask = (maskbg == 0);

