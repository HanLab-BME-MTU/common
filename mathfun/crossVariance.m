function [ Bcc ] = crossVariance(ts1,ts2,tFluc)
%function [ Bcc ] = crossVariance(sd,sCurForce_sd,tFluc) calculates
%cross-variance of the two time series, ts1 and ts2, with a time span of
%tFluc. 
%   input:      ts1:    1xN time series 1
%               ts2:    1xN time sereis 2
%               tFluc:  time span (default 11)
%   output:     Bcc:    1xN cross variance
%                       the first value will appear after tFluc/2
% Sangyoon Han 2018 June

if nargin<3
    tFluc = 11;
end
halfTfluc = floor(tFluc/2);
lastFrameCC = find(~isnan(ts1),1,'last') -halfTfluc;
firstFrameCC = find(~isnan(ts1),1) +halfTfluc;
Bcc = NaN(size(ts1));

for jj=firstFrameCC:lastFrameCC
    ts1_segment = ts1(jj-halfTfluc:jj+halfTfluc); avgTS1 = mean(ts1_segment);
    ts2_seg = ts2(jj-halfTfluc:jj+halfTfluc); avgTS2 = mean(ts2_seg);
    sigma2cc=1/tFluc*sum((ts1_segment-avgTS1).*(ts2_seg-avgTS2));
    Bcc(jj) = sigma2cc/sqrt(avgTS1*avgTS2);
end

end

