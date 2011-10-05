function [status,pValue]=permTest(pop1,pop2,alpha,varargin)
% permTest performs a permutation test for means, returning the one-sided p-value
%
% the idea is to sample n1 and n2 values (with replacement) from the union
% of the populations pop1 and pop2.  if the null hypothesis is true (that
% the one distributions are the same), then the mean(sample1)-mean(sample2)
% should cluster around zero after many iterations.
%
% the one-sided p-value is then the proportion of sampled permutations
% where the difference in means was greater than or equal to the absolute
% value of the difference between the population means.

% input check
ip = inputParser;
ip.addRequired('pop1',@isnumeric);
ip.addRequired('pop2',@isnumeric);
ip.addRequired('alpha',@isscalar);
ip.addOptional('tail','both',@(x) any(strcmpi(x,{'both','right'})));
ip.parse(pop1,pop2,alpha,varargin{:})
tail=ip.Results.tail;

nReps=1000;

n1=length(pop1);
n2=length(pop2);

% get union of the two populations
bigPop=[pop1(:); pop2(:)];

% get absolute value of the difference between the actual population means
deltaPop=abs(mean(pop1)-mean(pop2));

% sample both populations with replacement nReps times to get a
% distribution of the differences
delta=zeros(nReps,1);
for i=1:nReps
    s1=randsample(bigPop,n1,'true');
    s2=randsample(bigPop,n2,'true');

    delta(i)=mean(s1)-mean(s2);
end

% two-sided p-value: proportion of abs(delta) values greater than deltaPop
%pValue = sum(abs(delta)>deltaPop)/nReps;

% calculate the one-sided p-value
if strcmpi(tail,'both')
    pValue = 2*(1-normcdf(deltaPop,mean(delta),std(delta)));
else
    pValue = 1-normcdf(deltaPop,mean(delta),std(delta));
end

status = pValue<=alpha;
%figure; hist(delta,25)