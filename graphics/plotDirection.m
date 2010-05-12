function plotDirection(linesIn,varargin)
%PLOTCONTOURS plots the input 2D lines with an arrow indicating the direction of the curve
%
% plotDirections(linesIn)
% plotContours(linesIn,plotStyle1,plotStyle2,...)
%
% Description:
%   
%   Plots the input 2D lines and shows the direction the points occur
%   within the line by placing a small arrow at the first point.  
% 
% Input:
% 
%   linesIn - A single 2xM matrix containing the line to plot, or a cell
%             array of 2xM matrices if multiple lines are to be plotted.
% 
%   plotStyleString - Optional. A string (or strings) specifying the
%                     style/color to use when plotting the lines (same as
%                     with the plot command)
%
% Hunter Elliott
% 4/2010
% NOTE: Replaces my old plotContours.m function.
%

if ~iscell(linesIn)
    linesIn = {linesIn};
end

plotArgs = varargin;

hold on

%Plot the lines
cellfun(@(x)(plot(x(1,:),x(2,:),plotArgs{:})),linesIn)

%Plot the first point with arrow
cellfun(@(x)(quiver(x(1,1),x(2,1),x(1,2)-x(1,1),x(2,2)-x(2,1),5,plotArgs{:},'MaxHeadSize',10)),linesIn);
