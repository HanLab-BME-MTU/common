function  [outTS,newX] = gapInterpolation(TS,nPoint)

%This function interpolates NaN in a TS gap. It does not interpolate borders.
%USAGE:
%       workTS = gapInterpolation(TS,nPoint)
%
%Input:
%       TS     - Time Series (row Vector)
%       nPoint - maximum size of the gaps to be interpolated
%
%Output:
%       workTS - Time Series with closed gaps (row vector)
%       newX   - new x-axis vector 
%
% 
% This function does not interpolate points at the borders 
% Marco Vilela, 2012

 

nObs   = numel(TS);
xi     = find(isnan(TS));
outTS  = TS(:)';
workTS = TS(:)';
newX   = 1:nObs;

 
if (~isempty(xi) && numel(xi) < nObs) 

    

    nanB         = findBlock(xi,1);
    %If NaN blocks in the beginning or end of the TS, Delete it
    leftBorder   = nanB{1}(1)==1;
    rightBorder  = find(nanB{end}(end)==nObs);
    workTS([nanB{leftBorder};nanB{rightBorder*end}]) = [];

    %Indexes for real points

    numIdx       = find(~isnan(workTS));
    %Blocks of real points
    numB         = findBlock(numIdx,1);
    %Fusing blocks that will have an interpolated number in between
    fusingB      = find(cell2mat(cellfun(@(x,y)  y(1) - x(end) <= nPoint+1,numB(1:end-1),numB(2:end),'Unif',0)));


    if ~isempty(fusingB)

        fusedB       = findBlock(fusingB,1);
        %Fused blocks with gaps <= nPoint
        fusedPoint   = cellfun(@(x) cat(1,numB{[x;x(end)+1]}),fusedB,'Unif',0);
        %New x-axis
        newX         = (nanB{1}(1)==1)*nanB{1}(end) + cell2mat(cellfun(@(x) x(1):x(end),fusedPoint,'Unif',0));
        interpF      = @(x,y) interp1(x,y(x),x(1):x(end));
        %Interpolated points
        interpPoint  = cellfun(@(x) interpF(x,workTS),fusedPoint,'Unif',0);

        outTS(newX)  = cell2mat(interpPoint)';

    end

end