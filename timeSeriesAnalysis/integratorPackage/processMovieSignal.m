function processMovieSignal(movieObject,varargin)
% PROCESSMOVIESIGNAL process time series from the sampled maps
%
% SYNOPSIS processMovieSignal(movieObject,paramsIn)
%
% INPUT   
%   movieObject - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson, Jan 2012

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieObject', @(x) isa(x,'MovieObject'));
ip.addOptional('paramsIn',[], @isstruct);
ip.addParamValue('waitbar',[], @ishandle);
ip.parse(movieObject,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous  signal preprocessing process                                                                  
iProc = movieObject.getProcessIndex('SignalProcessingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieObject.processes_)+1;
    movieObject.addProcess(SignalProcessingProcess(movieObject,...
        movieObject.outputDirectory_));                                                                                                 
end

signalProc = movieObject.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(signalProc,paramsIn);

stack=dbstack;
if ~any(strcmp('Process.run',{stack(:).name}));
    signalProc.run();
    return;
end

%% --------------- Initialization ---------------%%
% Call all the tools execution function
for i=1:numel(p.processingTools)
    p.processingTools.function(movieObject,i,wtBarArgs{:});
end