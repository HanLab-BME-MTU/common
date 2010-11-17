% Francois Aguet, 16 November 2010

function lm = locmax3d(img, mask)

if numel(mask)==1
    wx = mask;
    wy = mask;
    wz = mask;
    if mod(wx,2)==0 || mod(wy,2)==0 || mod(wz,0)==0
        error('Mask dimensions must be odd integers');
    end
    mask = ones(wy,wx);
elseif numel(mask)==3
    wx = mask(1);
    wy = mask(2);
    wz = mask(3);
    if mod(wx,2)==0 || mod(wy,2)==0 || mod(wz,0)==0
        error('Mask dimensions must be odd integers');
    end    
    mask = ones(wy,wx);
end

[ny,nx,nz] = size(img);
lm2D = zeros(ny,nx,nz);

for z = 1:nz
    lm2D(:,:,z) = ordfilt2(img(:,:,z), wx*wy, mask);  
end

lm = zeros(ny,nx,nz);

b = (wz-1)/2;
for z = 1+b:nz-b
    lm(:,:,z) = max(lm2D(:,:,z-b:z+b), [], 3);
end

lm(lm~=img) = 0;

% also set xy-borders to zero
b = (wx-1)/2;
lm(:,[1:1+b end-b:end],:) = 0;
b = (wy-1)/2;
lm([1:1+b end-b:end],:,:) = 0;