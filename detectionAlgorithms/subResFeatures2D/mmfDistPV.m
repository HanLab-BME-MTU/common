function pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree)

pValue = zeros(numMaxima);

for k = 1 : numMaxima-1
    for j = k+1 : numMaxima
        
        %calculate distance between the 2 maxima
        x1_x2 = maximaPos(j,1) - maximaPos(k,1);
        y1_y2 = maximaPos(j,2) - maximaPos(k,2);
        distance = sqrt(x1_x2^2+y1_y2^2);
        
        %get the standard deviation in the distance
        j1 = 3*(j-1)+1;
        k1 = 3*(k-1)+1;
        stdDist = x1_x2^2*(varCovMat(j1,j1) + ...
            varCovMat(k1,k1) - 2*varCovMat(j1,k1)) ...
            + y1_y2^2*(varCovMat(j1+1,j1+1) + ...
            varCovMat(k1+1,k1+1) - 2*varCovMat(j1+1,k1+1)) ...
            + 2*x1_x2*y1_y2*(varCovMat(j1,j1+1) - ...
            varCovMat(j1,k1+1) - varCovMat(j1+1,k1) + ...
            varCovMat(k1,k1+1));
        stdDist = sqrt(stdDist)/distance;
        
        %1-sided t-test: H0: T=0, H1: T>0
        %calculate test statistic (t-distributed)
        testStat = distance/stdDist;
        
        %get p-value
        pValue(j,k) = 1-tcdf(testStat,numDegFree);
        
    end
end

%% ~~~ the end ~~~