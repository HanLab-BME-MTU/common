% [h coeff R xLR yLR]=plotLR(x,y,x_label,y_label)
% 
% Linear Regression plot to fit the line y = mx + b to data. It models the
% relationship between a dependent variable Y and an independent variables X
% Input:
% x -> dependent variable
% y -> independent variable 
% x_label -> x label
% y_label -> y label
% Output:
% h -> figure handle
% coeff -> coefficients [m b]
% R -> linear or correlation
% xLR -> x values of fitted line
% yLR -> y values of fitted line

function [h coeff R xLR yLR]=plotLR(x,y,x_label,y_label)

[coeff] = polyfit(x,y,1);
R = corr(x,y);

figure,
h = gcf;
hold on;
grid off;
set(gcf, 'Color','white');
set(gca,'FontSize',16,'LineWidth',3);
set(gca,'FontWeight','Bold','XColor','k','YColor','k','Color','w');
minval = min(min(x),min(y));
maxval = max(max(x),max(y));
axis ([minval maxval minval maxval])
set(gca,'FontSize',16)
plot(x,y,'ko','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');
xLR=(round(minval):round(maxval));
yLR=coeff(1).*xLR+ coeff(2);
plot(xLR,yLR,'-k','LineWidth',2);
xlabel(x_label);
ylabel(y_label);
temp = sprintf('y = %.2f x + %.2f',coeff(1),coeff(2));
text((minval + ((maxval-minval)/8)),(maxval - (maxval-minval)/8),temp,'FontSize',16,'Color','k')
temp = sprintf('R^2 = %.4f',R^2);
text((minval + ((maxval-minval)/5)),(maxval - (maxval-minval)/5),temp,'FontSize',16,'Color','k')
