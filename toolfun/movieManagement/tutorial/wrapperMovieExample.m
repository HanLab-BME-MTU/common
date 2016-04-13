function wrapperMovieExample(movieDataOrProcess)

%% Registration

% Get MovieData object and Process
[movieData, process] = getOwnerAndProcess(movieDataOrProcess,'ExampleProcess',true);

%% Input check

%% Input/output
p = process.getParameters();
outputFile = fullfile(p.OutputDirectory, 'output.mat');
if ~isdir(p.OutputDirectory), mkdir(p.OutputDirectory), end

% logging input
nChan = numel(movieData.channels_);
inFilePaths = cell(nChan, 1);
inFilePaths{1} = movieData.channels_(1).channelPath_;
process.setInFilePaths(inFilePaths);

% logging output
nChan = numel(movieData.channels_);
outFilePaths = cell(nChan, 1);
outFilePaths{1} = outputFile;
process.setOutFilePaths(outFilePaths);


%% Algorithm

output = zeros(movieData.nFrames_, 1);
for t = 1 : movieData.nFrames_
    I = movieData.getChannel(1).loadImage(t);
    
    % Algorithm
    output(t) = myFunction(I);
        
end
save(outputFile, 'output');

end