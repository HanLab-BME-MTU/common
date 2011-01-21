function [x,y] = boundaryGeometry(bs,s)
%BOUNDARYGEOMETRY generic geometry M file for use with PDE toolbox.
%
% [x,y] = boundaryGeometry(bs,s)
%
% This function is designed for use with solvePDEVectorBoundary.m and the
% PDE toolbox.
% 
%
% Hunter Elliott
% Re-written 8/2010
%

%We need the boundary coordinates to be a global variable, because the PDE
%toolbox won't pass any extra arguments and won't accept a function handle
%as an input
global OBJ_BOUND;

%Arbitrary? Pretty sure this is only for the piecewise case.
nBoundSeg = 10;

%If no arguments supplied, return the number of boundary segments. This is
%standard for geometry m file
if nargin == 0
    x = nBoundSeg;
    return
end

%Initialze d matrix which defines the boundary segments
d = nan(4,nBoundSeg);

%Check if curve runs clockwise or counterclockwise
tstCurve = ppval(OBJ_BOUND,linspace(0,1,5e2));
isClockwise = isCurveClockwise(tstCurve);

%Adjust labelling so that inside is always labelled 1
if isClockwise
    dirVec = [0 1];
else
    dirVec = [1 0];
end

for j = 1:nBoundSeg
    %The parameter runs from 0 to 1 along the cell border.
    %The interior is labelled 1 and the exterior zero
    d(:,j) = [ (j-1)/nBoundSeg j/nBoundSeg dirVec ]';        
end

%Standard output if 1 argument for geometry m-file
if nargin == 1
    x = d(:,bs);
    return
end

%Adjust so points are evenly spaced by arc length. Only necessary 
%if there are gaps in cell edge, remove if speed needed
iPars = linspace(0,1,OBJ_BOUND.pieces+1);
currEdge = fnval(OBJ_BOUND,iPars);
adjustedPar = pdearcl(iPars,currEdge,s,0,1);

% Evaluate this spline at the requested points 
%yy = ppval(OBJ_BOUND,s);
yy = ppval(OBJ_BOUND,adjustedPar);

x = nan(size(s));
y = nan(size(s));

%Split up the x & y coord for return
x(:) = yy(1,:);
y(:) = yy(2,:);


