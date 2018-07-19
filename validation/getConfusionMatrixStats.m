function [Acurracy TruePositiveRate TrueNegativeRate FalsePositiveFraction FalseNegativeFraction] = getConfusionMatrixStats(ConfusionMatrix)
% Usage 
% getConfusionMatrixStats(ConfusionMatrix);
%
%reference is webpage : http://www2.cs.uregina.ca/~hamilton/courses/831/notes/confusion_matrix/confusion_matrix.html 
%
% ConfusionMatrix has the structure 
%
%   TN   FP
%   FN   TP
%           	        Predicted 	
% 		             Negative	Positive
% Actual 	Negative	TN	      FP
% 	        Positive 	FN	      TP
%
%


TN = ConfusionMatrix(1,1);
FN = ConfusionMatrix(2,1);
FP = ConfusionMatrix(1,2);
TP = ConfusionMatrix(2,2);


%The accuracy (AC) is the proportion of the total number of predictions that were correct. It is determined using the equation: 
Acurracy = 100*(TN + TP)/(TN + FN + FP + TP);

%The true positive rate (TP) is the proportion of positive cases that were correctly identified, as calculated using the equation
TruePositiveRate = 100*TP/(TP + FN);

%The true negative rate (TN) is defined as the proportion of negatives cases that were classified correctly, as calculated using the equation: 
TrueNegativeRate = 100*TN/(TN + FP);


%The false positive rate (FP) is the proportion of negatives cases that were incorrectly classified as positive, as calculated using the equation: 
%FalsePositiveFraction = 100*FP/(TP + FN);
FalsePositiveFraction = 100*FP/(TN + FP);

%The false negative rate (FN) is the proportion of positives cases that were incorrectly classified as negative, as calculated using the equation: 
FalseNegativeFraction = 100*FN/(FN + TP);

%precision (P) is the proportion of the predicted positive cases that were correct, as calculated using the equation: 
Precision = 100*FP/(FP+TP);



