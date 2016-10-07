function [ objectFig ] = ftexplore( I )
%ftexplore Explore localized Fourier Transform
%
% INPUT
% I - image to explore
% OUTPUT
% objectFig - figure with the image
%
% For a real image, I, the real component of the Fourier coefficients
% contains information relative to the pixel I(1,1). We can visualize the
% Fourier transform relative to any point by circularly shifting the image
% such that the point of interest is at I(1,1).
%
% Furthermore, we can focus on information local to the pixel of the image
% by pixelwise multiplication of the image by a Gaussian distribution
% centered on that pixel. This is equivalent to convolving frequency space
% by the corresponding Gaussian.
%
% This tool allows you to visualize the component of the Fourier
% coefficients relative to and local to an arbitrary pixel in the image.
% This is especially useful for filtering in frequency space since the
% response of the filter is the pointwise multiplication with 

I = double(I);
Iz = I - mean(I(:));

objectFig = figure;
him = imshow(I,[]);
I_sz = size(I);
h_pt = impoint(gca,[I_sz(2)/2,I_sz(1)/2]);
h_ellipse = imellipse(gca,[-20+I_sz(2)/2 -20+I_sz(1)/2 41 41]);
h_ellipse.setFixedAspectRatioMode(true);
h_ellipse.addNewPositionCallback(@posCallback);


delta = zeros(I_sz);
delta(1) = 1;

posCallback(getPosition(h_ellipse));


freqFig = figure;
hfreq = imshow(abs(fftshift(fft(Iz))),[]);
colorbar
cm = vertcat(bsxfun(@times,(1-(0:32)/32)',[0 1 1]),bsxfun(@times,(1:32)'/32,[1 1 0]));

while(isvalid(h_ellipse))  
    showLocalizedFreq(getPosition(h_ellipse));
    wait(h_ellipse);
end


    function posCallback(p)
        center = [p(2)+p(4)/2,p(1)+p(3)/2];
        center = round(center);
        gaussianWindow = imgaussfilt(circshift(delta,center),p(3)/3,'FilterSize',round(p(3))*4+1);
        cdata = imfuse(gaussianWindow.*I,I,'Scaling','independent');
        set(him,'CData',cdata);
        setPosition(h_pt,fliplr(center));
    end
    function showLocalizedFreq(p)
        center = [p(2)+p(4)/2,p(1)+p(3)/2];
        center = round(center);
        gaussianWindow = ifftshift(imgaussfilt(fftshift(delta),p(3)/3,'FilterSize',round(p(3))*4+1));
        Is = circshift(Iz,-center+1);
        Is = Is.*gaussianWindow;
%         Is = Is - mean(Is(:));
        If = real(fftshift(fft2(Is)));
        cmax = max(abs(If(:)));
        figure(freqFig);
        set(hfreq,'CData',If);
%         imshow(real(fftshift(fft2(Is))),[]);
        caxis(hfreq.Parent,cmax*[-1 1]);
        colormap(hfreq.Parent,cm);
    end

end

