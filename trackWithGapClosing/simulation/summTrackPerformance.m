function [fracLinksCorrect,fracLinksWrong,fracGapsCorrect,fracGapsWrong,...
    fracLinksFP,fracGapsFP] = summTrackPerformance(resultsSim)

n3 = size(resultsSim,3);

if n3 == 1 %simulation with no false positives
    
    nMiss = size(resultsSim,1);
    
    [fracLinksCorrect,fracLinksWrong,fracGapsCorrect,fracGapsWrong] = ...
        deal(NaN(nMiss,1));
    fracLinksFP = [];
    fracGapsFP = [];
    
    for i = 1 : nMiss
        tmp = vertcat(resultsSim(i,:).linkStats0);
        fracLinksCorrect(i) = sum(tmp(:,3))/sum(tmp(:,1));
        fracLinksWrong(i) = sum(tmp(:,4))/sum(tmp(:,1));
    end
    
    for i = 1 : nMiss
        tmp1 = vertcat(resultsSim(i,:).gapStats0);
        tmp2 = vertcat(resultsSim(i,:).gapInfo);
        fracGapsCorrect(i) = length(find(tmp1(:,2)==0))/size(tmp2,1);
        fracGapsWrong(i) = length(find(tmp1(:,2)==1))/size(tmp2,1);
    end
    
else %simulation with false positives
    
    nFP = size(resultsSim,1);
    nMiss = size(resultsSim,2);
    
    [fracLinksCorrect,fracLinksWrong,fracGapsCorrect,fracGapsWrong,...
        fracLinksFP,fracGapsFP] = deal(NaN(nMiss,nFP));
    
    for i = 1 : nFP
        for j = 1 : nMiss
            tmp = vertcat(resultsSim(i,j,:).linkStats0);
            fracLinksCorrect(j,i) = sum(tmp(:,3))/sum(tmp(:,1));
            fracLinksWrong(j,i) = sum(tmp(:,4))/sum(tmp(:,1));
            fracLinksFP(j,i) = sum(tmp(:,5))/sum(tmp(:,1));
        end
    end
    
    for i = 1 : nFP
        for j = 1 : nMiss
            tmp1 = vertcat(resultsSim(i,j,:).gapStats0);
            tmp2 = vertcat(resultsSim(i,j,:).gapInfo);
            fracGapsCorrect(j,i) = length(find(tmp1(:,2)==0))/size(tmp2,1);
            fracGapsWrong(j,i) = length(find(tmp1(:,2)==1))/size(tmp2,1);
            fracGapsFP(j,i) = length(find(tmp1(:,2)==2))/size(tmp2,1);
        end
    end
    
end
