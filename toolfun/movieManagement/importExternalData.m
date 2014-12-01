function importExternalData(movieData, processClass, varargin)
% importExternalData imports detection results into the dummy process
%
% Sebastien Besson Nov 2014

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addRequired('processClass',@ischar);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,processClass,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous dummy detection processes
iProc = movieData.getProcessIndex(processClass,1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    processConstr = str2func(processClass);
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(processConstr(movieData,...
        movieData.outputDirectory_));
end
dummyProc = movieData.getProcess(iProc);

%Parse input, store in parameter structure
p = parseProcessParams(dummyProc,paramsIn);

%% --------------- Initialization ---------------%%

% Get channel paths and initialize process paths and output dirs
nChan = numel(movieData.channels_);

%
channelIndex = find(~cellfun(@isempty, p.InputData));
channelIndex = channelIndex(:)';

% Setup the  input directories
inFilePaths = cell(1,nChan);
for i = channelIndex
    assert(exist(p.InputData{i}, 'file') == 2 || ...
        exist(p.InputData{i}, 'dir') == 7);
    inFilePaths{1, i} = p.InputData{i};
end
dummyProc.setInFilePaths(inFilePaths);

% Setup the output directories
outFilePaths = cell(1,nChan);
if ~isdir(p.OutputDirectory), mkdir(p.OutputDirectory); end
for  i = channelIndex
    [~, name, ext] = fileparts(p.InputData{i});
    if isempty(ext)
        copyfile(p.InputData{i}, fullfile(p.OutputDirectory, name));
        outFilePaths{1,i} = fullfile(p.OutputDirectory, name);
    else
        copyfile(p.InputData{i}, p.OutputDirectory);
        outFilePaths{1,i} = fullfile(p.OutputDirectory, [name ext]);
    end
    
end
dummyProc.setOutFilePaths(outFilePaths);

disp('Finished importing detection!')