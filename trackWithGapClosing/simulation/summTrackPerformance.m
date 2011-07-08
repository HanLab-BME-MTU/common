function [fracLinksCorrect,fracLinksWrong,fracGapsCorrect,fracGapsWrong] = ...
    summTrackPerformance(resultsSim)

[fracLinksCorrect,fracLinksWrong,fracGapsCorrect,fracGapsWrong] = deal(NaN(8,1));

for i = 1 : 8
    tmp = vertcat(resultsSim(i,:).linkStats0);
    fracLinksCorrect(i) = sum(tmp(:,3))/sum(tmp(:,1));
    fracLinksWrong(i) = sum(tmp(:,4))/sum(tmp(:,1));
end

for i = 2 : 8
    tmp1 = vertcat(resultsSim(i,:).gapStats0);
    tmp2 = vertcat(resultsSim(i,:).gapInfo);
    fracGapsCorrect(i) = length(find(tmp1(:,2)==0))/size(tmp2,1);
    fracGapsWrong(i) = length(find(tmp1(:,2)==1))/size(tmp2,1);
end

