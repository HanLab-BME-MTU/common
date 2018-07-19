% [h,m,d,BAstats] = plotBA(var1,var2,x_label,y_label)
% 
% Bland-Altman plot to determine the agreement between two variables.
% Input:
% var1 -> variable 1
% var2 -> variable 2
% x_label -> x axis label
% y_label -> y axis label
% CI -> confidence interval
% Output:
% h -> figure handle
% m -> mean of two variables
% d -> difference between two variables
% BAstats -> structure that includes: 
%             *n = sample size
%             *bias = mean of the differences
%             *std = std of the differences
%             *CI = confidence interval for bias
%             *SE = standard error
%             *tstat = t-value
%             *DF = degrees of freedom
%             *p = p-value
%             *lowerCI = confidence interval for lower interval of agreement
%             *upperCI = confidence interval for upper interval of agreement
%
% Reference
% STATISTICAL METHODS FOR ASSESSING AGREEMENT BETWEEN TWO METHODS OF CLINICAL MEASUREMENT
% J. Martin Bland, Douglas G. Altman


function [h,m,d,BAstats] = plotBA(var1,var2,x_label,y_label)

n=max(size(var1));
m = (var1 + var2)/2;       % mean of two measurements
d = var1-var2;             % difference between two measurements
md = mean(d);         % mean of the differences (Bias)
sd = std(d);          % std of the differences

[nullh,pvalue,ci,dstats] = ttest(d); %One-sample t-test
tlmtagr=tinv(1-(0.05/2),dstats.df);
SE=sqrt((sd^2)/n);
SElmtagr=sqrt((3*(sd^2))/n);

md_01= ones(1,length(var1)).*md;
mp_2sd= md+ones(1,length(var1)).*2*sd;
mm_2sd= md-ones(1,length(var2)).*2*sd;

bias_cip = ones(1,length(var1))*ci(1);
bias_cim = ones(1,length(var2))*ci(2);
lower_cip = mm_2sd+ones(1,length(var1))*SElmtagr*tlmtagr;
lower_cim = mm_2sd-ones(1,length(var2))*SElmtagr*tlmtagr;
upper_cip = mp_2sd+ones(1,length(var1))*SElmtagr*tlmtagr;
upper_cim = mp_2sd-ones(1,length(var2))*SElmtagr*tlmtagr;

BAstats.n=n;
BAstats.bias=md;
BAstats.std=sd;
BAstats.CI=ci;
BAstats.SE=SE;
BAstats.tstat=dstats.tstat;
BAstats.df=dstats.df;
BAstats.p=pvalue;
BAstats.lowerCI=[lower_cim(1) lower_cip(1)];
BAstats.upperCI=[upper_cim(1) upper_cip(1)];

figure;
h = gcf;
box on;
set(gcf,'Color','w');
set(gca,'LineWidth',3,'FontSize',16);
set(gca,'FontWeight','Bold','XColor','k','YColor','k','Color','w');

hold on;
plot(m, md_01, '-k','LineWidth',3);     % bias
plot(m, mp_2sd, '-k','LineWidth',3);   % lower limit of agreement
plot(m, mm_2sd, '-k','LineWidth',3);   % upper limit of agreement
plot(m, d,'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',4);

plot(m, bias_cip, '-k','LineWidth',2);   
plot(m, bias_cim, '-k','LineWidth',2);
plot(m, lower_cip, '-k','LineWidth',2);   
plot(m, lower_cim, '-k','LineWidth',2);
plot(m, upper_cip, '-k','LineWidth',2);   
plot(m, upper_cim, '-k','LineWidth',2);

% % Linear regression on diff over mean
% [diff_fresult,diff_gof,diff_output] = fit(m',d','poly1');
% diff_p = confint(diff_fresult);
% liney=(diff_fresult.p1).*points + (diff_fresult.p2);
% plot(points,liney,'k-.','LineWidth',3);


space=(min(m)+max(m))/10;
axis ([min(m)-0.5*space max(m)+1.5*space -max(abs(d)) max(abs(d))])

temp = sprintf('%.2f',md);

text(max(m)+space/10,md,temp,'FontSize',16,'Color','k');
text(max(m)+space/10,(md+2*sd),'2 STD','FontSize',16,'Color','k');
text(max(m)+space/10,(md-2*sd),'-2 STD','FontSize',16,'Color','k');
    
xlabel(x_label,'FontSize',16,'FontWeight','bold','VerticalAlignment','cap');
ylabel(y_label,'FontSize',16,'FontWeight','bold','VerticalAlignment','bottom');
