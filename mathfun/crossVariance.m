function [ Bcc, bgBcc ] = crossVariance(ts1,ts2,tFluc)
%function [ Bcc ] = crossVariance(sd,sCurForce_sd,tFluc) calculates
%cross-variance of the two time series, ts1 and ts2, with a time span of
%tFluc. 
%   input:      ts1:    1xN time series 1
%               ts2:    1xN time sereis 2
%               tFluc:  time span (default 11)
%   output:     Bcc:    1xN cross variance
%                       the first value will appear after tFluc/2
%               bgBcc:  distribution of cross variance out of the 
% Sangyoon Han 2018 June

if nargin<3
    tFluc = 11;
end
halfTfluc = floor(tFluc/2);
lastFrameCC = min(find(~isnan(ts1),1,'last'),find(~isnan(ts2),1,'last')) -halfTfluc;
firstFrameCC = max(find(~isnan(ts1),1),find(~isnan(ts2),1)) +halfTfluc;
Bcc = NaN(size(ts1));
bgBccSeries = NaN(size(ts1));

pp=0;
% fluc1=zeros(1,(lastFrameCC-firstFrameCC+1)*tFluc);
% fluc2=zeros(1,(lastFrameCC-firstFrameCC+1)*tFluc);

for jj=firstFrameCC:lastFrameCC
    ts1_segment = ts1(jj-halfTfluc:jj+halfTfluc); avgTS1 = mean(ts1_segment);
    ts2_seg = ts2(jj-halfTfluc:jj+halfTfluc); avgTS2 = mean(ts2_seg);
    sigma2cc=1/tFluc*sum((ts1_segment-avgTS1).*(ts2_seg-avgTS2));
    Bcc(jj) = sigma2cc/sqrt(abs(avgTS1*avgTS2));
    
    % For background,
    % 1. get the average fluctuation from the two signals
    pp=pp+1;
    fluc1 = ts1_segment-avgTS1;
    fluc2 = ts2_seg-avgTS2;
%     fluc1(tFluc*(pp-1)+1:tFluc*(pp)) = ts1_segment-avgTS1;
%     fluc2(tFluc*(pp-1)+1:tFluc*(pp)) = ts2_seg-avgTS2;
    % 2. Make random time series out of distribution of fluc1 and fluc2
    randTS1 = std(fluc1)*randn(size(fluc1))+avgTS1;
    randTS2 = std(fluc2)*randn(size(fluc2))+avgTS2;
    % 3. Calculate the cross variance again out of randTS1 and randTS2
    ts1_segRand = randTS1; avgTS1Rand = mean(ts1_segRand);
    ts2_segRand = randTS2; avgTS2Rand = mean(ts2_segRand);
    sigma2ccRand=1/tFluc*sum((ts1_segRand-avgTS1Rand).*(ts2_segRand-avgTS2Rand));
    bgBccSeries(jj) = sigma2ccRand/sqrt(abs(avgTS1Rand*avgTS2Rand));
end

bgBcc = nanstd(bgBccSeries); %nanmean(bgBccSeries)+

% % 2. Make random time series out of distribution of fluc1 and fluc2
% randTS1 = std(fluc1)*randn([1 (lastFrameCC-firstFrameCC+1)*tFluc]);
% randTS2 = std(fluc2)*randn([1 (lastFrameCC-firstFrameCC+1)*tFluc]);

% % 3. Calculate the cross variance again out of randTS1 and randTS2
% bgBccSeries = NaN(1,pp);
% 
% for jj=1:pp
%     ts1_segment = randTS1(tFluc*(jj-1)+1:tFluc*(jj)); avgTS1 = mean(ts1_segment);
%     ts2_seg = randTS2(tFluc*(jj-1)+1:tFluc*(jj)); avgTS2 = mean(ts2_seg);
%     sigma2cc=1/tFluc*sum((ts1_segment-avgTS1).*(ts2_seg-avgTS2));
%     bgBccSeries(jj) = sigma2cc/sqrt(abs(avgTS1*avgTS2));
% end
end

