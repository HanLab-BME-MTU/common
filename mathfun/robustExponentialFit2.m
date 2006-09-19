function [u, sigmaU, goodIdx, plotAx] = robustExponentialFit2(Y, A, verbose)
%ROBUSTEXPONENTIALFIT2 robustly fits exponential functions in a more stable way than robusExponentialFit, but without additive constant
%
% SYNOPSIS: [u, sigmaU, goodIdx] = robustExponentialFit2(Y, A)
%
% INPUT Y: n-by-1 vector of function values. The code cannot handle mixed
%           positive and negative values of Y. In that case, it will call
%           robustExponentialFit, and there is no possiblility of
%           complicated A's.
%		A: (opt) design matrix. Default is [ones(size(Y)), (1:n)']. For
%           multiple processes with the same exponential decay, A has
%           multiple columns 1:end-1, with A(i,j)=1 where Y(i) is
%           a result from process j.
%           Please always use the last column of A for the time/x values,
%           as the code might otherwise return erroneous results!
%       verbose: (opt) Plots the exponential fit. Default: 0. Can also be a
%           handle to the plot axes.
%
% OUTPUT u: parameters of the exponential. With the default A, the
%           exponential form will be y=u(1)*exp(u(2)*t).
%        sigmaU: uncertainty in u
%		 goodIdx: List of inlier rows in Y
%        plotAx: if verbose, handle to the plot axes. [] otherwise.
%
% REMARKS: The function will perform a robust linear fit to log(Y) with
%          weights Y.
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 03-Jul-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=====================
%% TEST INPUT
%=====================

% check for verbose
def_verbose = 0;
if nargin < 3 || isempty(verbose)
    verbose = def_verbose;
end

% Check for correct size of Y
if nargin == 0 || isempty(Y)
    error('robustExponentialFit needs at least a non-empty vector of function values as input!')
end

Y = returnRightVector(Y);



% Check for A
lengthY = length(Y);
if nargin < 2 || isempty(A)
    % default A
    A = [ones(lengthY,1), (1:lengthY)'];
else
    A = returnRightVector(A,lengthY,'r');
end

% check for bad rows
badYRows = ~isfinite(Y) | Y==0;
Y(badYRows) = [];
A(badYRows,:) = [];

% reset lengthY
lengthY = length(Y);

if isempty(Y)
    error('all entries of Y are either 0, NaN or Inf!')
end

% check for mixed positive and negative values
useREF = false;
if all(Y<0)
    % all Y negative. Simply divide by -1
    ySign = -1;
    Y = -Y;
elseif any(Y<0) && any(Y>0)
    % use robustExponentialFit
    useREF = true;
    ySign = 1;

    % check for complicated A
    if size(A,2) > 2
        error('cannot combine mixed sign Y-values with multiprocess A')
    end
else
    % all is well
    ySign = 1;
end

%===========================




%===========================
%% FIT
%===========================

if useREF
    % use robustExponentialFit
    [u, sigmaU, goodIdx] = robustExponentialFit(A(:,end),Y,0,0);

else
    % transform data and do weighted linear fit

    B = log(Y);
    W = Y;
    [u, sigmaU, goodIdx] = linearLeastMedianSquares(A,B,1./W);

    % transform u, sigmaU of multiplicative constants (! do proper error
    % propagation!)
    u(1:end-1) = exp(u(1:end-1));
    sigmaU(1:end-1) = u(1:end-1) .* sigmaU(1:end-1);

end


% transform sign if necessary
u(1:end-1) = ySign * u(1:end-1);


%================================



%===============================
%% PLOT
%==============================

if verbose
    
    if floor(verbose)~= verbose && ishandle(verbose)
        % axes handle has been given
        plotAx = verbose;
    else
        % create new figure, new axes
        figure('Name','Exponential Fit');
        plotAx = gca;
    end
    hold on
    
    
    t = A(:,end);
    colorOrder = get(gca,'ColorOrder');
    
    yFit=exp(A * [log(u(1:end-1));u(end)]);
    for i = 1:size(A,2)-1
        pIdx = find(A(:,i));
        % make sure that the time is ordered when plotting the estimate
        plotData = [t(pIdx),yFit(pIdx)];
        plotData = sortrows(plotData,1);
        plot(plotData(:,1),plotData(:,2),'Color',colorOrder(wraparound(i,[1;size(colorOrder,1)]),:))        
        plot(t(pIdx),Y(pIdx),'.','Color',colorOrder(wraparound(i,[1;size(colorOrder,1)]),:))

    end
    badIdx = setdiff(1:lengthY,goodIdx);
    plot(t(badIdx), Y(badIdx),'*r');
    
else
    plotAx = [];
end