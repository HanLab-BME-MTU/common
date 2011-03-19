function [cleanContours,iClean] = cleanUpContours(contoursIn,nPtsMin)
%CLEANUPCONTOURS removes redundant points and short contours from the input contours
%
% cleanContours = cleanUpContours(contoursIn)
% cleanContours = cleanUpContours(contoursIn,nPtsMin)
% [cleanContours,iClean] = cleanUpContours(contoursIn,nPtsMin)
% 
% Description:
% 
%   The contours returned by contours.m or contourc.m should first be
%   seperated into a cell array using separateContours.m. These contours
%   will often contain "redundant" points - multiple points in the same
%   pixel of the image which was contoured. (I'm not sure why this is...)
%   This function removes these points. Also, they may contain very short
%   contours which surround a single pixel. These contours will be removed
%   if they are shorter than nPtsMin. Finally, contours which intersect
%   themselves will be seperated into two, non-intersecting contours. This
%   happens only in rare circumstances, where the isovalue passes directly
%   through a saddle point.
% 
% 
% Input:
% 
%   contoursIn - A cell array containing the contours to process, as output
%                by separateContours.m
% 
%   nPtsMin - The minimum number of points a contour must contain to be
%             kept. Contours shorter than this will be removed. 
%             Optional. Default is 3 - just enough to make a triangle. 
% 
% 
% Output:
% 
%   cleanContours - The contours with redundant points and short contours
%                   removied.
% 
%   iClean - The indices of the returned clean contours in the original
%            contour array.
% 
% 
% 
% 
% Hunter Elliott 
% Re-Written 4/2010

%% ------- Input ------- &&

if nargin < 1 || isempty(contoursIn) || ~iscell(contoursIn)
    error('1st input "contoursIn" must be a cell-array of contours!')    
end

if nargin < 2 || isempty(nPtsMin)
    nPtsMin = 3; %Just enough to make a triangle - only contours that can actually contain something will be returned.
end

%% ----- Parameters ----- %%

%threshold for considering two points on a contour redundant
distThreshold = .5;

nContours = length(contoursIn);

%Minimum length of contours to consider splitting
nMinSplit = nPtsMin*3;

%% ------ Clean Up----- %%

%Find the sum difference in the coordinates at each point
%Don't use distance, to keep it fast!
allDiffs = cellfun(@(x)(sum(vertcat(abs(diff(x(1,:))),abs(diff(x(2,:)))))),contoursIn,'UniformOutput',false) ;

%Remove points which are seperated by less than the threshold
cleanContours = arrayfun(@(x)(contoursIn{x}(:,[(allDiffs{x} > distThreshold) true])),1:nContours,'UniformOutput',false)';

%Find the contours that are long enough
nPall = cellfun(@(x)(size(x,2)),cleanContours);
iClean = find(nPall >= nPtsMin);

%Return only these long contours
cleanContours = cleanContours(iClean)';
nPall = nPall(iClean);


%Check if the contour intersects itself
i1 = cell(1,nContours);
i2 = cell(1,nContours);
%Only check contours larger than nMinSplit
[~,~,i1(nPall>nMinSplit),i2(nPall>nMinSplit)]  = cellfun(@(x)(...
                              intersections(x(1,:),x(2,:))),...
                              cleanContours(nPall>nMinSplit),'UniformOutput',false);

%Exclude the first point because closed contours meet there.
i1 = cellfun(@(x)(x(x~=1)),i1,'UniformOutput',false);
i2 = cellfun(@(x)(x(x~=1)),i2,'UniformOutput',false);
                          
iSelfInt = find(cellfun(@(x)(~isempty(x)),i1));

if ~isempty(iSelfInt)

    iClean = arrayfun(@(x)(x),iClean,'UniformOutput',false)';
    
    for j = iSelfInt

        c1 = cleanContours{j}(:,max(ceil(i1{j})):min(floor(i2{j})));
        c2 = cleanContours{j}(:,[1:min(floor(i1{j})) max(ceil(i2{j})):end]);

        %This is a little trick to keep the indices correct/contours in
        %order
        cleanContours{j} = cell(1,2);
        cleanContours{j}{1} = c1;
        cleanContours{j}{2} = c2;
        iClean{j} = ones(1,2) .* iClean{j};        
    end
    cleanContours = cat(2,cleanContours{:});
    iClean = cat(2,iClean{:});
end

%Get rid of contours that were split into 'too small' contours
nPall = cellfun(@(x)(size(x,2)),cleanContours);
cleanContours = cleanContours(nPall>=nPtsMin)';
iClean = iClean(nPall>=nPtsMin);
