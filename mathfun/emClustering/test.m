% firstfilename=('D:\s400\actin\S400ACTIN001.tif');
% 
% 
% [fileName,dirName] = uigetfile('*.tif','Select image')
% firstfilename=[dirName,fileName];

LOWER_CUT=0.005;
UPPER_CUT=0.99;
   
   
   
img=imread('D:\s400\actin\S400ACTIN001.tif');
img_hist=imhist(img,16383);


%take xx% percent of the pixels
total_pix=sum(img_hist);
%find 3% 
i=1;
while sum(img_hist(1:i)) < LOWER_CUT*total_pix
   i=i+1;
end
i_3=i;
%find 97%
i=1;
while sum(img_hist(1:i)) < UPPER_CUT*total_pix
   i=i+1;
end
i_97=i;

img_hist(1:i_3)=0;
img_hist(i_97:length(img_hist))=0;


[ii jj]=find(img_hist);
img_hist(ii(length(ii)):length(img_hist))=[];
img_hist(1:ii(1))=[];
t0=round(length(img_hist)/2);

%generate points
ind=1;
for i=1:6:length(img_hist)
    for j=1:4000:img_hist(i)
        A(ind,1)=i;
        A(ind,2)=j;
        ind=ind+1;
    end
end
figure,plot(A(:,1),A(:,2),'.');


