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
info = @real;
hax(1) = viscircles(hfreq.Parent,repmat([101 101],4,1),[10 20 40 80],'DrawBackgroundCircle',false,'EdgeColor',[1 1 1],'LineStyle',':','LineWidth',1);
hax(2) = line(xlim-0.5,repmat(sum(xlim)/2,1,2),'Color',[1 1 1],'LineWidth',1,'LineStyle',':','Parent',hfreq.Parent);
hax(3) = line(repmat(sum(ylim)/2,1,2),ylim-0.5,'Color',[1 1 1],'LineWidth',1,'LineStyle',':','Parent',hfreq.Parent);
menu = uicontextmenu(freqFig);
m_real = uimenu(menu,'Label','Real','Checked','on','Callback',@selectInfo);
m_imag = uimenu(menu,'Label','Imag','Checked','off','Callback',@selectInfo);
m_abs  = uimenu(menu,'Label','Magnitude','Checked','off','Callback',@selectInfo);
m_angle = uimenu(menu,'Label','Phase','Checked','off','Callback',@selectInfo);
m_axes = uimenu(menu,'Label','Axes','Checked','on','Callback',@toggleAxes,'Separator','on');
hfreq.UIContextMenu = menu;
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
        If = info(fftshift(fft2(Is)));
        cmax = max(abs(If(:)));
        figure(freqFig);
        set(hfreq,'CData',If);
%         imshow(real(fftshift(fft2(Is))),[]);
        caxis(hfreq.Parent,cmax*[-1 1]);
        colormap(hfreq.Parent,cm);
    end
    function selectInfo(hObject,~)
        set(menu.Children,'Checked','off');
        switch(hObject.Label)
            case 'Real'
                info = @real;
                set(m_real,'Checked','on');
            case 'Imag'
                info = @imag;
                set(m_imag,'Checked','on');
            case 'Magnitude'
                info = @abs;
                set(m_abs,'Checked','on');
            case 'Phase'
                info = @angle;
                set(m_angle,'Checked','on');
        end
        showLocalizedFreq(getPosition(h_ellipse));
    end
    function toggleAxes(hObject,~)
        switch(hax(1).Visible)
            case 'on'
                set(hax,'Visible','off');
                m_axes.Checked = 'off';
            case 'off'
                set(hax,'Visible','on');
                maxes.Checked = 'on';
        end
    end

end

