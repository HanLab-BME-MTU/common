function sem = stdErrMean(xvector)
% function sem = stdErrMean(xvector) calculates standard error of mean of
% xvector.
% Sangyoon Han 2016 March
sem = nanstd(xvector)./sqrt(sum(~isnan(xvector)));