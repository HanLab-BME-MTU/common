function[cpar,pvr,dpvr, cpar2]=clusterQuantRipley(mpm,imsizex,imsizey);
% clusterQuantRipley calculates a quantitative clustering parameter based on
% Ripley's K-function (a spatial statistics function for point patterns)
%
% SYNOPSIS   [cpar,pvr,dpvr]=clusterQuantRipley(mpm,imsizex,imsizey);
%       
% INPUT      mpm:   mpm file containing (x,y) coordinates of points in the
%                   image in succesive columns for different time points
%            imsizex:   x-size of the image (maximum possible value for x-coordinate)
%            imsizey:   y-size of the image (maximum possible value for
%                       y-coordinate)
%            cf:  correction =1 circumferebnce +2 area
%            NOTE: in Johan's mpm-files, the image size is 1344 x 1024
%               pixels
%            NOTE2: this function uses Ripley's circumference correction
%            NOTE3: Although this function is not actually named after 
%                   Lt. Ellen Ripley, she certainly would deserve to have 
%                   a kick-ass matlab function named after her.
%
%
% OUTPUT     cpar:  for each time point, a single clustering parameter 
%                   value is extracted from the pvr function  
%            pvr:   for each plane (time point), the function calculates
%                   the function pvr=points vs radius, i.e. the number of
%                   points contained in a circle of increasing radius
%                   around an object, averaged over all objects in the
%                   image
%                   NOTE: the default size for the radius implicit in the 
%                   pvr function is [1,2,3...,minimsize] where minimsize is the
%                   smalller dimension of imsizex,imsizey
%                   This function corresponds to Ripley's K-function (/pi)
%            dpvr:  H(r) function
%            cpar2: additional parameters
%
% DEPENDENCES   ClusterQuantRipley uses {pointsincircle,clusterpara}
%               ClusterQuantRipley is used by { }
%
% Dinah Loerke, October 7th, 2004

%create vector containing x- and y-image size
matsiz=[imsizex imsizey];
rs=round(min(matsiz)/2);

%determine size of mpm-file
[nx,ny]=size(mpm);

%initialize results matrix pvr; x-dimension equals the employed number of
%values for the circle radius, y-dimension equals number of planes of the
%input mpm-file
pvr=zeros(rs,(ny/2));
dpvr=zeros(rs,(ny/2));
cpar=[1:(ny/2)];
cpar2=[1:(ny/2)];


%initialize temporary coordinate matrix matt, which contains the object 
%coordinates for one plane of the mpm
matt=zeros(nx,2);

%cycle over all planes of series, using two consecutive columns of mpm input
%matrix as (x,y) coordinates of all measured points
for k=1:(round(ny/2))
    
    %matt is set to two consecutive columns of input matrix m1
    matt(:,:)=mpm(:,(2*k-1):(2*k));
    
    %since the original mpm file contains a lot of zeros, these zeros are 
    %deleted in the temporary coordinate matrix to yield a matrix containing
    %only the nonzero points of matt, smatt
    nz1=size(nonzeros(matt(:,1)));
    nz2=size(nonzeros(matt(:,2)));
    if(nz1==nz2)
        smatt=[nonzeros(matt(:,1)), nonzeros(matt(:,2)) ];
    else
        disp(['different number of point entries for x and y in plane ',num2str(k)]);
    end
    
    %comment/uncomment the next five lines if you want to monitor progress
    %prints number of objects for every 10th line
    if(mod(k,10)==0)
        [smx,smy]=size(smatt);
        tempnp=max([smx,smy]);
        disp('  plane  number of objects');
        tempi=[k, tempnp];
        disp(tempi);
    end
    
    %now determine number of objects in circle of increasing radius,
    %averaged over all objects in smatt, and normalized with point density
    %tempnp/(msx*msy)
    [pvrt]=pointsincircle(smatt,matsiz);
    %result is already normalized with point density tempnp/(msx*msy)
    pvr(:,k)=pvrt(:);
    
    %from the calculated function pvrt (number of points versus circle
    %radius), calculate a quantitative clustering parameter, cpar
    %somewhat arbitrarily defined as positive integral
    [cpar(k),dpvr(:,k), cpar2(k)]=clusterpara(pvrt);
        
end
%normalize cpar with initial value
ini=(cpar(1)+cpar(2))/2;
cpar=cpar/ini;
end


function[cpar1, dpvrt, cpar2]=clusterpara(pvrt);
%clusterpara calculates a quantitative cluster parameter from the input
%function (points in circle) vs (circle radius)
% SYNOPSIS   [cpar]=clusterpara(pvrt);
%       
% INPUT      pvrt:   function containing normalized point density in 
%                   circle around object
%                   spacing of points implicitly assumes radii of 1,2,3...
%   
% OUTPUT     cpar:    cluster parameter
%            dpvrt:   difference function of p vs r
%
% DEPENDENCES   clusterpara uses {DiffFuncParas}
%               clusterpara is used by {FractClusterQuant}
%
% Dinah Loerke, September 13th, 2004

%calculate difference L(d)-d function, using L(d)=sqrt(K(d))
%since K(d) is already divided by pi
len=max(size(pvrt));
de=(1:len);
%diff=sqrt(abs(pvrt))-de;
% H(r) funtcion from poission clustering
Hr=pvrt-de.^2;
dpvrt=Hr;
%extract parameters from diff
[cpar1,cpar2]=DiffFuncParas(Hr);

end    

    
function[p1,p2]=DiffFuncParas(Hr);
%DiffFuncParas calculates a number of quantitative cluster parameter 
%from the input function, the difference function
% SYNOPSIS   DiffFuncParas(diff);
%       
% INPUT      diff:  difference function as calculated in clusterpara
%                   vector with len number of points
%%   
% OUTPUT     p1,p2:  cluster parameters
%                    currently: p1=inclination of the first rise of the
%                                   H(t) function (clustering) 
%                               of total clustering)
%                               p2= position of first rise
%
% DEPENDENCES   DiffFuncParas uses {}
%               DiffFuncParas is used by {clusterpara}
%
% Dinah Loerke, September 13th, 2004

%% firstpoint: point where diff function systematically rises above zero 
%% definition: diff>0 AND (diff)'>0 to exlude noisy one-point rises above
%% zero
vec=Hr;
dvec=diff(vec);
detvec=vec;
detvec(vec<0)=0;
detvec(dvec<0)=0;
    
firstpoint=min(find(detvec));

% beginning point begp and end point endp for calculating inclination of the
% function diff
begp=round(firstpoint+3);
endp=begp+10;
%% these values can be refined to include variable lengths according to end of first
%% steep rise (considering point where inbclination starts to drop)
%% this will be done in next version

inclination=sum(dvec(begp:endp))/(endp-begp);
p1=inclination;
p2=firstpoint;

% OLD VERSION
% len=max(size(diff));
% p1=0;
% p2=0;
% for i=1:len
%     if(diff(i)>0)
%         p1=p1+diff(i);
%         if(diff(i)==max(diff))
%             p2=diff(i);
%         end
%     end
% end
end



function[m2]=pointsincircle(m1,ms)
%pointsincircle calculates the average number of points in a circle around
%a given point as a function of the circle radius (averaged over all points
%and normalized by total point density); this function is called Ripley's
%K-function in statistics, and is an indication of the amount of clustering
%in the point distribution
% 
% SYNOPSIS   [m2]=pointsincircle(m1,ms);
%       
% INPUT      m1:   matrix of size (n x 2) containing the (x,y)-coordinates of n
%                  points
%            ms: vector containing the parameters [imsizex imsizey] (the 
%                   x-size and y-size of the image)
%            NOTE: in Johan's mpm-files, the image size is 1344 x 1024
%               pixels
%
%
% OUTPUT     m2:    vector containing the number of points in a circle 
%                   around each point, for an increasing radius;
%                   radius default values are 1,2,3,....,min(ms)
%                   function is averaged over all objects in the
%
% DEPENDENCES   pointsincircle uses {distanceMatrix, circumferenceCorrectionFactor}
%                   (distanceMatrix, circumferenceCorrectionFactor added to this file)
%               pointsincircle is used by {FractClusterQuant }
%
% Dinah Loerke, October 4th, 2004


[lm,wm]=size(m1);

%for points at the edges (where the circle of increasing size is cut off by
%the edges of the image), this function corrects for the reduced size of 
%the circle using the 
%function circumferenceCorrectionFactor

msx=ms(1);
msy=ms(2);
minms=min(ms);
rs=round(minms/2);

%create neighbour matrix m3
%matrix m3 contains the distance of all points in m1 from all points
%in itself
[mdist]=distanceMatrix(m1,m1);

%create numpoints vector (number of points in circle of corresponding radius)
%loop over all radius values between 1 and minms
%initialize m2 vector
%m2=1:rs;

%create corrections factor matrix (same dimension as mdist)
%contains correction factor for precise radii (point distances) around each
%point; for the zero entry at identity (p11,p22,p33), cfm equals one
corrFacMat=ones(lm);
for n=1:lm
    corrFacMat(n,:) = circumferenceCorrectionFactor(m1(n,1),m1(n,2),mdist(n,:),msx,msy);
end
thresh_mdist = mdist;

for r=rs:-1:1
    %for given radius, set all values of mdist higher than the radius value
    %to zero
    
    thresh_mdist(thresh_mdist > r) = 0;
    
    %count all leftover points equally => set to one
    %to zero
    thresh_mdistones = thresh_mdist;
    thresh_mdistones(thresh_mdist > 0) = 1;
    
    %weigt every counted point with the circumference correction factor
    %calculated previously in corrFacMat
    tempfinal=thresh_mdistones./corrFacMat;    
    %sum over entire matrix to get number of points
    npv=sum(tempfinal(:));
        
    %to average, divide sum by number of points (=columns)
    npv=(npv/lm);
    
    %in order to be able to quantitatively compare the clustering in 
    %distributions of different point densities, this npv value must now 
    %be corrected for overall point density, which is lm/msx*msy; the 
    %resulting normalized function is (if we also divide by pi to scale for
    %the circle area) more or less a simple square function;
    %it is a perfect square function for a perfectly random distribution of
    %points
    m2(r)=npv/(pi*(lm-1)/(msx*msy));
    %using (lm-1) and not lm is Marcon&Puech's correction (2003)
end

end
  

function[m2]=distanceMatrix(c1,c2)
%this subfunction makes a neighbour-distance matrix for input matrix m1
%input: c1 (n1 x 2 points) and c2 (n2 x 2 points) matrices containing 
%the x,y coordinates of n1 or n2 points
%output: m2 (n1 x n2) matrix containing the distances of each point in c1 
%from each point in c2
[ncx1,ncy1]=size(c1);
[ncx2,ncy2]=size(c2);
m2=zeros(ncx1,ncx2);
for k=1:ncx1
    for n=1:ncx2
        d=sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);
        m2(k,n)=d;
    end
end
end
    
function[corfac]=circumferenceCorrectionFactor(xx,yy,rr,msx,msy)
%circumference correction calculates a vector containing the correction factor
%(for edge correction in Ripley's k-function) for values of rr
%circumference correction: fraction of circumference of circle centered at
%point P=(xx,yy) with radius rr (inside rectangular image) falling into the 
%rectangle - this fraction becomes smaller as the point gets closer to one
%of the rectangle's edges, and as the radius of the circle increases
%if the circle falls completely inside the rectangle, the value is zero

%1. this function assumes that rr is a vector

% SYNOPSIS   [corfac]=circumferenceCorrectionFactor2(xx,yy,rr,msx,msy)
%       
% Dinah Loerke, October 6, 2004


x=min(xx,(msx-xx));
y=min(yy,(msy-yy));

rmax=max(size(rr));
corfac=ones(rmax,1);

for i=1:rmax
    r=rr(i);
    %if both x and y are smaller than r
    if(r>0)
        
        if((x<r)&&(y<r))
            if( (x<r) && (y<r) && (sqrt(x^2+y^2)>r) )
                corfac(i)=(2*asin(x/r)+2*asin(y/r))/(2*pi);
            else
                corfac(i)=(0.5*pi+asin(x/r)+asin(y/r))/(2*pi);
            end
         %if either x OR y OR neither is smaller than r 
        else
             z=min( min(x,y),r );
            corfac(i)=(pi+2*asin(z./r))/(2*pi);
        end
    end
end

end

