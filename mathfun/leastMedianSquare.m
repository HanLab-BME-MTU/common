function [u,sigma0,ssq]=leastMedianSquare(functionString,u0,options,parameters)
%leastMedianSquare calculates the least median squares for a function handed down as string using parameters (including x/y-data) from the structure parameters
%
%SYNOPSIS [x,sigma]=leastMedianSquare(functionString,x0,options,parameters)
%
%INPUT functionString       string which specifies the function to be minimized, e.g. '(y-(u(1)*x+u(2)))^2'
%      u0                   vector of the initial values for the unknowns in the function
%      options              option structure for optimization routines. Set them with optimset
%      parameters           structure containing all parameters of the function e.g. y and x
%
%AS TESTING THE INPUT FOR SPELLING ERRORS IN THE FUNCTION STRING OR THE FIELDNAMES OF PARAMETERS IS FAR TOO COMPLICATED
%TO DO WITHIN THE FUNCTION, PLEASE BE CAREFUL WITH YOUR INPUT!
%
%OUTPUT  u                  vector with estimates for the unknowns u, e.g. a and b
%        sigma              estimate for sigma0 of the data according to the fit
%        ssq                sum of squares (y-y°)^2
%
%REMEMBER TO CHECK THE NUMBER OF LOCAL MINIMA IN THE REGION OF YOUR SOLUTION
%
%note that you can do a weighted fitting specifying the function as
%(y-f(x))*P*(y-f(x))', where P is the weight matrix

%c: Jonas, 11/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
parmString=[];
k=3; %value important for calculation of sigma, see Danuser, 1992 or Rousseeuw & Leroy, 1987
magicNumber2=1.4826^2; %see same publications

%get names of parameters for fminsearch
parameterNames=fieldnames(parameters);
parameterNamesN=length(parameterNames);


%prepare fminsearch
funStr=['''median(',functionString,')'', ''u'''];
funStrShort=['median(',functionString,')'];

for i=1:parameterNamesN
    parmString=[parmString,', ',parameterNames{i}];
    eval([parameterNames{i},'=parameters.',parameterNames{i},';']);
    funStr=[funStr,', ''',parameterNames{i},''''];
end

eval(['u=fminsearch(inline(',funStr,'),u0,options ',parmString,');']);

%calculate statistics
eval(['medRes2=', funStrShort,';']);
eval(['res2=', functionString,';']);

%testvalue to calculate weights
testValue=res2/(magicNumber2*medRes2);

%goodRows: weight 1, badRows: weight 0
goodRows=find(testValue<=k^2);

ssq=sum(res2);
sigma0=sqrt(sum(res2(goodRows))/(length(goodRows)-4));

