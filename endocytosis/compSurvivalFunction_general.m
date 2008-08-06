function [survFunc,pvec] = compSurvivalFunction_general(data1, field1, data2, field2)
% compSurvivalFunction_general compares the survival functions for
% different conditions
%
% SYNOPSIS [survFunc,pvec] = compSurvivalFunction_general(data1, field1,
% data2, field2)
%
% INPUT     data1    = experiment structure with first data set
%           field1  = field in data1 from which the lifetime histogram is
%                   read, e.g. 'survivalFunction_InRegion'
%           data2   = experiment structure with second data set
%           field2  = field in data1 from which the lifetime histogram is
%                   read, E.G. 'survivalFunction_OutRegion'
%
% OUTPUT    survFunc = survival function 
%                       (first row DATA1, second row Data2)
%           pvec    = vector of p-values;
%                     1: KS-test on difference between data1 and data2
%                     2: 5% threshold for KS-test on bootstrapped
%                     subsamples of data1
%                     3. 5% threshold for KS-test on bootstrapped
%                     subsamples of data2 
%
% NOTE: The current version of the function assumes that ALL MOVIES in the
% data structure are acquired at the same framerate, or that the user
% chooses to treat them as being the same, and that the survival functions
% are normalized to the value for the lifetime 1 frame.
%
%
% last modified DATE: 25-Jun-2008 (Dinah)
% last modified DATE: 23-Jul-2008 (Dinah)
% last modified DATE: 06-Aug-2008 (Dinah)


% averaging can only be performed up until the minimum common length of all 
% movies, so we determine the shortest movie length in the structure
n1 = length(data1);
n2 = length(data2);

for i=1:n1
    mlvec1(i) = data1(i).movieLength;
end
minlen1 = min(mlvec1);

for i=1:n2
    mlvec2(i) = data2(i).movieLength;
end
minlen2 = min(mlvec2);

minlen = min(minlen1,minlen2);


histmat1 = nan*zeros(length(data1),minlen);
histmat2 = nan*zeros(length(data2), minlen);

% read data for each movie 
for k=1:n1
    
    if isfield(data1,field1)
        survFun = getfield(data1(k), field1);
        % normalize individual survival functions to minlen length
        survFunNorm = survFun(1:minlen)/sum(survFun(1:minlen));
        histmat1(k,:) = survFunNorm;
    else
        error(['function requires a structure field called ',field1]);
    end          
end

% read data for each movie 
for k=1:n2
    
    if isfield(data2,field2)
        survFun = getfield(data2(k), field2);
        survFunNorm = survFun(1:minlen)/sum(survFun(1:minlen));
        histmat2(k,:) = survFunNorm;
    else
        error(['function requires a structure field called ',field2]);
    end          
end



% average survival functions
survFun1_AVE = nanmean(histmat1,1);
survFun2_AVE = nanmean(histmat2,1);

% define output
survFunc(1,:) = survFun1_AVE;
survFunc(2,:) = survFun2_AVE;

% compare averages
[H,pval_av] = kstest2(survFun1_AVE,survFun2_AVE);


%% bootstrap

% number of bootstrap runs
nbs = 1000;

% bootstrap first data set
for b=1:nbs
    % position vector
    pos = randsample(n1,n1,true); 
    % local bootstrap set from normalized matrix
    cmat = histmat1(pos,:);    
    % average
    bootstrapAVE = nanmean(cmat,1); 
    bootstrap1mat(b,:) = bootstrapAVE;   
    % comp between average and subsampled average
    [H,pval_bs] = kstest2(bootstrapAVE,survFun1_AVE);
    bootstrap1_pval(b) = pval_bs;

end

for b=1:nbs
    % position vector
    pos = randsample(n2,n2,true); 
    % local bootstrap set from normalized matrix
    cmat = histmat2(pos,:);    
    % average
    bootstrapAVE = nanmean(cmat,1); 
    bootstrap2mat(b,:) = bootstrapAVE;
    % comp between average and subsampled average
    [H,pval_bs] = kstest2(bootstrapAVE,survFun2_AVE);
    bootstrap2_pval(b) = pval_bs;
end

% bootstrap confidence intervals
sort_bs1 = sort(bootstrap1_pval);
sort_bs2 = sort(bootstrap2_pval);

cflevel1 = sort_bs1(round(0.05*nbs));
cflevel2 = sort_bs2(round(0.05*nbs));


% output
pvec(1) = pval_av;
pvec(2) = cflevel1;
pvec(3) = cflevel2;


%% display
figure; hold on;
plot(survFun1_AVE,'r-'); 
plot(survFun2_AVE,'b-');

xlabel('lifetime (frames)');
ylabel('survival function (fraction)');
legend('data1','data2');
axismax = 1.1*max(max(survFun1_AVE),max(survFun2_AVE));
axis([0 minlen -0.001 axismax]);

textstr1 = ['KS-test avdiff 1-2 p=',num2str(pval_av)];
textstr2 = ['KS-test 5% boostrap1 p=',num2str(cflevel1)];
textstr3 = ['KS-test 5% boostrap2 p=',num2str(cflevel2)];

text(50,0.9*axismax,textstr1);
text(50,0.85*axismax,textstr2);
text(50,0.8*axismax,textstr3);

end % of function







    
