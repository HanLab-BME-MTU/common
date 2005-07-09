function[acnx, acny]=normACfunc(image)
%[acnx, acny]=normACfunc(image)
%creates normalized spatial autocorrelation function (of mean zero) of
%image with wraparound, one in x-, one in y-direction
%INPUT: image = should be double precision; if image is type uint8 or uint16,
%       use normACfunc(double(image)) instead
%OUTPUT: acnx = normalized autocorrelation function of image projcted in x
%               -direction
%        acny = normalized autocorrelation function of image projcted in y
%               -direction
%last modified by Dinah 07/08/2005

[sizex,sizey]=size(image);
dbimage=double(image);
dbimageCon=dbimage';
im_xproject=dbimage(:);
im_yproject=dbimageCon(:);

%initialize
shiftedvecx=im_xproject;
shiftedvecy=im_yproject;
acnx=im_xproject;
acny=im_yproject;
%the mean m ist subtracted from the image further down to ensure that the
%mean of the autocorrelation is zero
m=mean(im_xproject);

%variable for progress monitoring
p=max(size(im_xproject));
ct=100*round(p/1000);

h = waitbar(0,'This might take awhile...');
%loop over r = spatial shift in pixels
for r=0:(p-1)
    waitbar(r/(p-1))
    %shiftedvec equals original vector shifted by r, with wraparound
    shiftedvecx(1:(p-r))=im_xproject((1+r):p);
    shiftedvecx((p-r):p)=im_xproject(1:(1+r));
    
    shiftedvecy(1:(p-r))=im_yproject((1+r):p);
    shiftedvecy((p-r):p)=im_yproject(1:(1+r));
    
    %acn =  normalized autocorrelation function int(f(s)*f(s+r))/int(f(s)^2)
    %the total mean of the vector is subtracted to ensure that the
    %autocorrelation function has mean zero (and thus drops to zero for
    %perfectly uncorrelated images)
    acnx(r+1)=sum((im_xproject-m).*(shiftedvecx-m))/sum((im_xproject-m).^2);
    acny(r+1)=sum((im_yproject-m).*(shiftedvecy-m))/sum((im_yproject-m).^2);

end % of for

close(h);
account=acnx;
account=0:(p-1);
plot(account,acnx,'r-');
my=min(min(acnx),min(acny));
mx=4*max(sizex,sizey);
axis([-10 mx my 1]);
hold on
plot(account,acny,'b-');
hold off

