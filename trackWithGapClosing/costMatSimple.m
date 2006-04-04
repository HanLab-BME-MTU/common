function [costMat,noLinkCost,errFlag] = costMatSimple(movieInfo,costMatParams)
%COSTMATSIMPLE provides a simple cost matrix for linking features between 2 time points
%
%SYNOPSIS function [costMat,noLinkCost,errFlag] = costMatSimple(movieInfo,costMatParams)
%
%INPUT  movieInfo    : A 2x1 array (corresponding to the 2 time points of 
%                      interest) containing the fields:
%             .xCoord      : Image coordinate system x-coordinate of detected
%                            features [x dx] (in pixels).
%             .yCoord      : Image coorsinate system y-coordinate of detected
%                            features [y dy] (in pixels).
%             .amp         : Amplitudes of PSFs fitting detected features [a da].
%       costMatParams: Structure with the following fields:
%             .searchRadius: Maximum distance between two features in two
%                            consecutive time points that allows linking 
%                            them (in pixels).
%             .maxAmpRatio : Maximum ratio between the amplitudes of two
%                            features in two censecutive time points that 
%                            allows linking them.
%             .noLnkPrctl  : Percentile used to calculate the cost of
%                            linking a feature to nothing. Use -1 if you do
%                            not want to calculate this cost.
%
%OUTPUT costMat      : Cost matrix.
%       noLinkCost   : Cost of linking a feature to nothing, as derived
%                      from the distribution of costs.
%       errFlag      : 0 if function executes normally, 1 otherwise.
%
%REMARKS The cost for linking feature i in time point t to feature j 
%in time point t+1 is given by
%(distance(i,j)^2)x(max(amp(i),amp(j))/min(amp(i),amp(j))).
%
%Khuloud Jaqaman, March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

costMat = [];
noLinkCost = [];
errFlag = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatSimple')
    disp('--costMatSimple: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

searchRadius = costMatParams.searchRadius;
maxAmpRatio = costMatParams.maxAmpRatio;
noLnkPrctl = costMatParams.noLnkPrctl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cost matrix calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of features in the 2 time points
n = movieInfo(1).num;
m = movieInfo(2).num;

%replicate x,y-coordinates at the 2 time points to get n-by-m matrices
x1 = repmat(movieInfo(1).xCoord(:,1),1,m);
y1 = repmat(movieInfo(1).yCoord(:,1),1,m);
x2 = repmat(movieInfo(2).xCoord(:,1)',n,1);
y2 = repmat(movieInfo(2).yCoord(:,1)',n,1);

%calculate the square distances between features in time points t and t+1
costMat = (x1-x2).^2 + (y1-y2).^2;

%assign NaN to all pairs that are separated by a distance > searchRadius
indx = find(costMat>searchRadius^2);
costMat(indx) = NaN;

%replicate the feature amplitudes at the 2 time points to get n-by-m
%matrices
a1 = repmat(movieInfo(1).amp(:,1),1,m);
a2 = repmat(movieInfo(2).amp(:,1)',n,1);

%divide the larger of the two amplitudes by the smaller value
ampRatio = a1./a2;
for j=1:m
    for i=1:n
        if ampRatio(i,j) < 1
            ampRatio(i,j) = 1/ampRatio(i,j);
        end
    end
end

%assign NaN to all pairs whose amplitude ratio is larger than the
%maximum allowed
indx = find(ampRatio>maxAmpRatio);
ampRatio(indx) = NaN;

%multiply the distance between pairs with the ratio between their
%amplitudes
costMat = costMat.*ampRatio;

%determine noLinkCost
if noLnkPrctl ~= -1
    noLinkCost = prctile(costMat(:),noLnkPrctl);
end

%replace NaN, indicating pairs that cannot be linked, with -1
costMat(find(isnan(costMat))) = -1;


%%%%% ~~ the end ~~ %%%%%

