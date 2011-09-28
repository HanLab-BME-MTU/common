function movieData = getMovieMasksMSS(movieData,varargin)
%THRESHOLDMOVIE applies multi-scale steerable filtering to every frame in input movie
%
% movieData = thresholdMovie(movieData,paramsIn)
%
% Applies manual or automatic thresholding to every frame of the input
% movie and then writes the resulting mask to file as a binary .tif in a
% sub-folder of the movie's analysis directory named "masks"
%
% Input:
% 
%   movieData - A MovieData object describing the movie to be processed, as
%   created by setupMovieDataGUI.m
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
% 
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
% Output:
%
%   movieData - the updated MovieData object with the thresholding
%   parameters, paths etc. stored in it, in the field movieData.processes_.
%
%   The masks are written to the directory specified by the parameter
%   OuptuDirectory, with each channel in a separate sub-directory. They
%   will be stored as binary, bit-packed, .tif files. 
%
%
% Sebastien Besson, Sep 2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous threshold processes from this function                                                                              
iProc = movieData.getProcessIndex('MSSSegmentationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(MSSSegmentationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end

segProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(segProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows'),
    wtBar = waitbar(0,'Initializing...','Name',segProc.getName());
else
end

if ~all(segProc.checkChanNum(p.ChannelIndex))
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

%Read various constants
imDirs  = movieData.getChannelPaths();
imageFileNames = movieData.getImageFileNames();
nFrames=movieData.nFrames_;

% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = imDirs{i};
end
segProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outputDir = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outputDir{1,i} = [p.OutputDirectory filesep 'mask_for_channel_' num2str(i)];
    mkClrDir(outputDir{1,i})
end
segProc.setOutFilePaths(outputDir);

%% ---------------Mask calculation ---------------%%% 


disp('Starting applying MSS segmentation...')
%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);

% Anonymous functions for reading input/output
inImage=@(chan,frame) [imDirs{chan} filesep imageFileNames{chan}{frame}];
outMask=@(chan,frame) [outputDir{chan} filesep 'mask_' numStr(frame) '.tif'];


logMsg = @(chan) ['Please wait, applying multi-scales filter to channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;
for i=1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    disp(logMsg(iChan))
    disp(imDirs{iChan});
    disp('Results will be saved under:')
    disp(outputDir{iChan});
    
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iChan)); end
    
    for j=1:nFrames
        % Read image apply mask and save the output
        currImage = double(imread(inImage(iChan,j)));
        mask=getCellMaskMSS(currImage,'Scales',p.Scales,'FilterOrder',p.FilterOrder);        
        imwrite(mask,outMask(iChan,j));

        % Update the waitbar
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
    
end

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished segmenting!')
