function [u,goodRows,sigma0]=leastMedianSquares(functionString,u0,options,parameters)
%leastMedianSquare calculates the least median squares for a function handed down as string using parameters (including x/y-data) from the structure parameters
%
%SYNOPSIS [u,goodRows,sigma0]=leastMedianSquares(functionString,u0,options,parameters)
%
%INPUT functionString       string which specifies the function to be
%                           minimized, e.g. '(y-(u(1)*x+u(2)))'. The unknown has to be 'u'. Do not square! 
%      u0                   vector of the initial values for the unknowns in the function
%      options              option structure for optimization routines. Set them with optimset
%      parameters           structure containing all parameters of the function e.g. y and x
%
%AS TESTING THE INPUT FOR SPELLING ERRORS IN THE FUNCTION STRING OR THE FIELDNAMES OF PARAMETERS IS FAR TOO COMPLICATED
%TO DO WITHIN THE FUNCTION, PLEASE BE CAREFUL WITH YOUR INPUT!
%
%OUTPUT  u                  vector with estimates for the unknowns u, e.g. a and b
%        goodRows           rows in the data that are inliers
%        sigma              estimate for std of error of the data according to the fit
%
%REMEMBER TO CHECK THE NUMBER OF LOCAL MINIMA IN THE REGION OF YOUR
%SOLUTION (fminsearch is sensitive to intial conditions)
%
% Possible improvements: - Use weights (e.g. field in parameter structure),
%                          and build the objective function anlalogous to
%                          linearLeastMedianSquares
%                       -  Allow boundaries (switch to fminbnd in that case)
%
%c: Jonas, 3/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
parmString=[];
k=3; %value important for calculation of sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magicNumber2=1.4826^2; %see same publications

u0 = u0(:);

%get names of parameters for fminsearch
parameterNames=fieldnames(parameters);
parameterNamesN=length(parameterNames);

%prepare fminsearch
funStringFmin = ['(' functionString,').*(',functionString ')']; % we want median(fs.*fs), not median fs'*fs!!
funStringMedian=['''median(',funStringFmin,')'', ''u'''];


for i=1:parameterNamesN
    parmString=[parmString,', ',parameterNames{i}]; % write variable into parameter string
    eval([parameterNames{i},'=parameters.',parameterNames{i},';']); % assign value to variable
    funStringMedian=[funStringMedian,', ''',parameterNames{i},'''']; % write variable into inline object
end

u=eval(['fminsearch(inline(',funStringMedian,'),u0,options ',parmString,');']);
uMedian = u;

%calculate statistics
eval(['res2=', funStringFmin,';']);
medRes2 = median(res2);

%testvalue to calculate weights
testValue=res2/(magicNumber2*medRes2);

%goodRows: weight 1, badRows: weight 0
goodRows=find(testValue<=k^2);

ssq=sum(res2);
sigma0=sqrt(sum(res2(goodRows))/(length(goodRows)-4));


% write function string for least squares formulation
funStringLsq = ['''(' functionString ')''''*(' functionString ')'', ''u'' ']; %now we want the sum of funStringFmin

% update parameters: take only good rows
for i=1:parameterNamesN
    eval([parameterNames{i},'=parameters.',parameterNames{i},'(goodRows);']); % reassign
    funStringLsq=[funStringLsq,', ''',parameterNames{i},'''']; % write variable into inline object
end



u=eval(['fminsearch(inline(',funStringLsq,'),u0,options ',parmString,');']);
    