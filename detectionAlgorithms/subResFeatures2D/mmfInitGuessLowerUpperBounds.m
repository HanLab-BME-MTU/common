function [x0,lb,ub] = mmfInitGuessLowerUpperBounds(maximaPosT,maximaAmpT,...
    bgAmpT,psfSigma,clusterPixels,firstFit)

%feature positions
x0 = maximaPosT; %initial guess
%         lb = repmat(min(clusterPixels),numMaximaT,1); %lower bound
%         ub = repmat(max(clusterPixels),numMaximaT,1); %upper bound
%         lb = x0-1; %lower bound
lb = x0 - 2*psfSigma; %lower bound
minPos = min(clusterPixels);
lb(lb(:,1)<minPos(1),1) = minPos(1);
lb(lb(:,2)<minPos(2),2) = minPos(2);
if ~firstFit
    lb(end,:) = minPos;
end
%         ub = x0+1; %upper bound
ub = x0 + 2*psfSigma; %upper bound
maxPos = max(clusterPixels);
ub(ub(:,1)>maxPos(1),1) = maxPos(1);
ub(ub(:,2)>maxPos(2),2) = maxPos(2);
if ~firstFit
    ub(end,:) = maxPos;
end
%feature amplitudes
x0 = [x0 maximaAmpT];
%         lb(:,3) = 1e-5;
%         ub(:,3) = 2*x0(:,3);
lb(:,3) = eps;
ub(:,3) = 1;
%background intensity
x0 = x0';
x0 = [x0(:); bgAmpT];
lb = lb';
%         lb = [lb(:); 1e-5];
lb = [lb(:); eps];
ub = ub';
%         ub = [ub(:); 2*bgAmpMax];
ub = [ub(:); 1];
