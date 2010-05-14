function varargout = parseProcessParams(procOb,paramList,paramIn)
%PARSEPROCESSPARAMS parses the parameters of the input process object
% 
% [param1,param2,...] = parseProcessParams(procObj)
% [param1,param2,...] = parseProcessParams(procObj,paramList)
% [param1,param2,...] = parseProcessParams(procObj,paramList,paramIn)
% 
% This function returns the individual parameters that are specified in the
% input process object's field "funParams" If the parameter list
% "paramList" is specified, then each of the parameters in this list will
% be returned, and if one is missing and exception will be thrown. If the
% input parameters "paramIn" is specified, these will override or
% substitute for any values in procOb.funParams_
% 
% Input:
%   
%   procOb - An object of the class Process.m
% 
%   paramList - 1xM or Mx1 cell array of character arrays specifying the names of
%   the M parameter fields which are expected to be present in
%   procOb.funParams_
% 
%   paramIn - 1xN or Nx1 cell array, where N is an even number, containing
%   parameter name/value pairs where element i is the parameter name and
%   element i+1 is the value for that parameter.
%
% Hunter Elliott
% 5/2010
%


if nargin < 1 || isempty(procOb) || ~isa(procOb,'Process')
    error('You must input a process object as the first input!')
end

if ~isstruct(procOb.funParams_) || ~isobject(procOb.funParams_) 
    error('The funParams_ field of the input object must be a structure or object containing fields with parameter values!')
end

if nargin < 2 || isempty(paramList)
    %If no parameter list was provided, just use every field in the
    %funParams
    paramList = fieldnames(procOb.funParams_);    
end

if nargin < 3 || isempty(paramIn)
    paramIn = {};
end

nPar = numel(paramList);

if ~isEven(nPar)
    error('The input paramers must be option name/value pairs!')
end

for j = 1:nPar

    %Check for this parameter in the input parameter list
    iPar = find(strcmpi(paramList{j},paramIn),1);
    
    %If the parameter was input, this overides existing parameters
    if ~isempty(iPar)
        procOb.funParams_.(paramList{j}) = paramList{j+1};
        varargout{j} = paramList{j+1};
    else
        %Look for this parameter in the object
        if isfield(procOb.funParams_,paramList{j})
            varargout{j} = procOb.funParams.(paramList{j});            
        else
            error(['Parameter ' paramList{j} ' was not input! Check process object and/or input parameters!'])
        end
    end        
end