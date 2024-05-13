function [insideIdx] = maskVectors(X,Y,myMask)
% function [insideIdx] = maskVectors(X,Y,myMask) finds indices in [X Y]
% points which are inside myMask.
% Example:
% [xmat ymat]=meshgrid(1:10, 1:10);
% myCenterPoint = [4,4];
% myR=3;
% myMask = sqrt((xmat-myCenterPoint(1)).^2+(ymat-myCenterPoint(2)).^2)<=myR;
% xVec = xmat(:); yVec = ymat(:);
% myIndexInside = maskVectors(xVec,yVec,myMask);
% imshow(myMask)
% hold on
% plot(xVec(myIndexInside),yVec(myIndexInside),'ro')
% Sangyoon Han
Npoints = length(X);
insideIdx = false(Npoints,1);

for ii=1:Npoints
    if round(Y(ii)) <= size(myMask,1) && ...
            round(Y(ii)) > 0 && ...
            round(X(ii)) <= size(myMask,2) && ...
            round(X(ii)) > 0 && ...
            myMask(round(Y(ii)),round(X(ii)))
        insideIdx(ii) = true;
    end
end
