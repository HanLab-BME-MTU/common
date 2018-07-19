function cMatrix  = getConfusionMatrix(imA,imB)
% Usage 
% cMatrix = getconfusionMatrix(predict,actual);
%
% CFmatrix 2std 
%
%   TN   FP
%   FN   TP
%           	        Predicted 	
% 		             Negative	Positive
% Actual 	Negative	TN	      FP
% 	        Positive 	FN	      TP
%
%

imA = logical(imA);
imB = logical(imB);


andAB = imA & imB;
% TruePositive TP 
%TP = sum(sum(andAB));
TP = sum(andAB(:));

orAB  = imA | imB;
% False Negative FN
%TN = sum(sum(~orAB));
TN = sum((~orAB(:)));

% False Positive
temp = orAB - imB;
FP = sum(temp(:));

% True Positive
temp = orAB - imA;
%FN = sum(sum(orAB - imA));
FN = sum(temp(:));

cMatrix = [TN FP ; FN TP];