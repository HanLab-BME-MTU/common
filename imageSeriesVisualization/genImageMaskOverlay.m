function [ imRGB ] = genImageMaskOverlay( im, mask, maskColor, maskAlpha )

    imr = im2uint8( mat2gray( im ) );
    img = imr;
    imb = imr;
    mask = mask > 0;
    maxVal = 255;
    
    imr(mask) = uint8( double( (1 - maskAlpha) * imr(mask) ) + maxVal * maskAlpha * maskColor(1) );
    img(mask) = uint8( double( (1 - maskAlpha) * img(mask) ) + maxVal * maskAlpha * maskColor(2) );
    imb(mask) = uint8( double( (1 - maskAlpha) * imb(mask) ) + maxVal * maskAlpha * maskColor(3) );
    
    imRGB = cat(3, imr, img, imb );
    
end