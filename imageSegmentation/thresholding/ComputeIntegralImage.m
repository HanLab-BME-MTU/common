function [imIntegral] = ComputeIntegralImage( im )

    imIntegral = im;
    for i = 1:ndims(im)
        imIntegral = cumsum(imIntegral,i);
    end

end