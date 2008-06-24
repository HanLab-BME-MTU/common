function likelihood = neg2LnLikelihoodXMEXcall(paramV,prob)


prob.user.xOrder = length(paramV) - (prob.user.arOrder + prob.user.maOrder) - 1;

sum1 = 0;
sum2 = 0;
numMissing = 0;
totalObs = 0;


for j = 1:length(prob.user.trajOut)

   
    trajLength = size(prob.user.trajOut(j).observations,1);
    
     if prob.user.xOrder < 1
        prob.user.trajIn(j).observations = zeros(trajLength,2);   
     end
    
    prob.user.TRAJ = cat(3, prob.user.trajIn(j).observations , prob.user.trajOut(j).observations);
    prob.user.TRAJ = permute( prob.user.TRAJ, [1 3 2]);
    prob.user.wnVariance = 1;
    probCall = prob.user;

    [sum1Tmp,sum2Tmp, numMissingTmp] = neg2LnLikelihoodXMEX(paramV,probCall);
    
    sum1 = sum1 + sum1Tmp;
    sum2 = sum2 + sum2Tmp;
    numMissing = numMissing + numMissingTmp;
    totalObs = totalObs + trajLength;    

end

likelihood = (sum1 + (totalObs - numMissing)*log(sum2));