function runDummySegmentation(movieData,varargin)
% runDummyDetection imports the detection results into the dummy process
%
% Sebastien Besson Nov 2014

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous dummy detection processes                                                                              
iProc = movieData.getProcessIndex('DummySegmentationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(DummySegmentationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
dummyProc = movieData.getProcess(iProc);

%Parse input, store in parameter structure
p = parseProcessParams(dummyProc,paramsIn);

%% --------------- Initialization ---------------%%

% Get channel paths and initialize process paths and output dirs
nChan = numel(movieData.channels_);

% Setup the  input directories
inFilePaths = cell(1,nChan);
for i = p.ChannelIndex;    
    assert(exist(p.InputData{i}, 'dir') == 7);
    inFilePaths{1, i} = p.InputData{i};
end
dummyProc.setInFilePaths(inFilePaths);

% Setup the output directories
outFilePaths = cell(2,nChan);
if ~isdir(p.OutputDirectory), mkdir(p.OutputDirectory); end
for  i = p.ChannelIndex;
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i)];
    copyfile([p.InputData{i} filesep '*'], outFilePaths{1,i});
end
dummyProc.setOutFilePaths(outFilePaths);

disp('Finished importing segmentation!')