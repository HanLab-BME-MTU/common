function paramOut = parseProcessParams(procOb,paramIn)
%PARSEPROCESSPARAMS parses the parameters of the input process object
% 
% [paramOut] = parseProcessParams(procObj,paramIn)
% 
% This function returns the individual parameters that are specified in the
% input process object's field "funParams" If the input parameter structure
% "paramIn" is specified, parameters specified in this structure will
% override any values in procOb.funParams_. They will be stored in the
% procOb and returned as outputs.
% 
% Input:
%   
%   procOb - An object of the class Process.m
% 
%   paramIn - A structure or object containing fields whose name's match
%   with the elements of paramList and who'se values are to be used as
%   parameters, replacing those stored in procOb. Optional. If not input,
%   all the values in procOb.funParams_ will be returned instead.
%
% Output:
%
%   paramOut - A structure containing the new parameters. Additionally the
%   parameters will be updated in the object.
%
% Hunter Elliott
% 5/2010
%


if nargin < 1 || isempty(procOb) || ~isa(procOb,'Process')
    error('You must input a process object as the first input!')
end

if nargin < 2
    paramIn = [];
end

if ~(isstruct(procOb.funParams_) || isobject(procOb.funParams_)) 
    error('The funParams_ field of the input object must be a structure or object containing fields with parameter values!')
end

%Get the list of parameter names from the process object
paramList = fieldnames(procOb.funParams_);    


nPar = numel(paramList);

%If the parameter was input, this overides existing parameters
paramOut = paramIn;

for j = 1:nPar
    %If the parameter wasn't input, use the default from the process.    
    if ~isfield(paramOut,paramList{j});                
        %Add this default value to the structure
        paramOut.(paramList{j}) = procOb.funParams_.(paramList{j});                    
    end        
end
procOb.setPara(paramOut);