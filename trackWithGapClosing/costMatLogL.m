function [costMat,errFlag] = costMatLogL(movieInfo,costMatParams)
%COSTMATLOGL provides a cost matrix for linking features between 2 time points based on amplitude and displacement statistics
%
%SYNOPSIS function [costMat,errFlag] = costMatLogL(movieInfo,costMatParams)
%
%INPUT  movieInfo    : A 2x1 array (corresponding to the 2 time points of 
%                      interest) containing the fields:
%             .xCoord      : Image coordinate system x-coordinate of detected
%                            features [x dx] (in pixels).
%             .yCoord      : Image coorsinate system y-coordinate of detected
%                            features [y dy] (in pixels).
%             .amp         : Amplitudes of PSFs fitting detected features [a da].
%       costMatParams: Structure with the following fields:
%             .dispSqLambda: Parameter of the exponential distribution
%                            that describes the displacement of a feature
%                            between two consecutive frames.
%             .ampDiffStd  : Standard deviation of the change in a feature's 
%                            amplitude between consecutive time points.
%             .maxDispSq   : Maximim squared displacement between consecutive 
%                            frames that allows the linking of two features.
%             .maxAmpDiffSq: Maximum squared change in amplitude between
%                            consecutive framces that allows the linking of
%                            two features.
%
%OUTPUT costMat      : Cost matrix.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The cost for linking feature i in time point t to feature j 
%in time point t+1 is given by -log[p(dI)p(dispSq)], where 
%dI = I(j;t+1) - I(i,t) is assumed to be normally distributed with mean 0
%and standard deviation ampDiffStd (supplied by user), and 
%dispSq = square of distance between feature i at t and feature j at t+1
%is assumed to be exponentially distributed with parameter dispSqLambda
%(supplied by user).
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costMat = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatLogL')
    disp('--costMatLogL: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

dispSqLambda = costMatParams.dispSqLambda;
ampDiffStd = costMatParams.ampDiffStd;
maxDispSq = costMatParams.maxDispSq;
maxAmpDiffSq = costMatParams.maxAmpDiffSq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of features in the 2 time points
n = movieInfo(1).num;
m = movieInfo(2).num;

%replicate the x,y-coordinates and amplitude at the 2 time points to get n-by-m matrices
x1 = repmat(movieInfo(1).xCoord(:,1),1,m);
y1 = repmat(movieInfo(1).yCoord(:,1),1,m);
a1 = repmat(movieInfo(1).amp(:,1),1,m);
x2 = repmat(movieInfo(2).xCoord(:,1)',n,1);
y2 = repmat(movieInfo(2).yCoord(:,1)',n,1);
a2 = repmat(movieInfo(2).amp(:,1)',n,1);

%calculate the squared displacement of features from t to t+1
dispSq = (x1-x2).^2 + (y1-y2).^2;
dispSq(find(dispSq>maxDispSq)) = NaN;

%calculate the squared change in the amplitude of features from t to t+1
ampDiffSq = (a1 - a2).^2;
ampDiffSq(find(ampDiffSq>maxAmpDiffSq)) = NaN;

%calculate the cost matrix
costMat = ampDiffSq/2/ampDiffStd^2 + dispSqLambda*dispSq;

%replace NaN, indicating pairs that cannot be linked, with -1
costMat(find(isnan(costMat))) = -1;


%%%%% ~~ the end ~~ %%%%%

