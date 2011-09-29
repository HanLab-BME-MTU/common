function distMat = bwMaxDirectDist(mask)

%Ill explain later... For now, just trust me.

[M,N] = size(mask);

distMat = zeros([size(mask) 4],'single');

distMat(1,:,1) = M;
distMat(M,:,2) = M;
distMat(:,1,3) = N;
distMat(:,N,4) = N;

for m = 2:M       
    distMat(m,~mask(m,:),1) = distMat(m-1,~mask(m,:),1) + 1;       
end
for m = M-1:-1:1   
    distMat(m,~mask(m,:),2) = distMat(m+1,~mask(m,:),2) + 1;       
end
for n = 2:N   
    distMat(~mask(:,n),n,3) = distMat(~mask(:,n),n-1,3) + 1;       
end
for n = N-1:-1:1   
    distMat(~mask(:,n),n,4) = distMat(~mask(:,n),n+1,4) + 1;       
end

%SET VALUES >= IMAGE DIMENSION EQUAL TO EACH OTHER!!
