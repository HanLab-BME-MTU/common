function sampleMovieSegmentation(movieData)
%SAMPLEMOVIESEGMENTEDAREAS samples the movie images in the segmented region(s)
% 
% sampleMovieSegmentedAreas(movieData)
% sampleMovieSegmentedAreas(movieData,paramsIn)
% 
% This function goes through each frame in the movie and samples the image
% intensities or process outputs in each segmented area. Various statistics
% for the pixels inside each segmented area are calculated, and written to
% disk. 
%
% **Note:** For a binary masks, one sample will be generated per iamge. To
% sample multiple areas independently an index image must be created, which
% each region to be sampled independently labelled witha  different integer
% value.
%
% 
% Input:
% 
%   movieData - A MovieData object describing the movie to be processed,
%   **AND** with a process object containing the parameters to use. See
%   SegmentationSamplingProcess.m for more info.
%
% Output:
%
%   Output samples will be written to disk in a sub-folder of the
%   MovieData's output directory, and all parameters logged in the
%   MovieData's process object.
% 

% Hunter Elliott
% 4/2013


%% ---------- Input ---------------- %%

if nargin < 1 || ~isa(movieData,'MovieData')   
    error('The first input must be a valid MovieData object!');        
end

   
%Check if the movie has a sampling process
iProc = movieData.getProcessIndex('SegmentationSamplingProcess',1,false);

assert(~isempty(iProc),'You must create the process object and add it to the MovieData before calling this function!')

sampProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(sampProc);

%Get the segmentation process index
if isempty(p.SegProcessIndex)    
    p.SegProcessIndex = movieData.getProcessIndex('MaskProcess',1,~p.BatchMode);                    
end
segProc = movieData.getProcess(p.SegProcessIndex);

if isempty(p.MaskChannelIndex)
    %Just use the first one that has valid masks
    p.MaskChannelIndex = find(segProc.checkChannelOutput,1,'first');
    assert(~isempty(p.MaskChannelIndex),'No channels have valid masks! Please segment the movie before attempting sampling!')
else    
    %Make sure the sementation has been successfully completed for the
    %specified channel
    assert(numel(p.MaskChannelIndex) == 1 && segProc.checkChannelOutput(p.MaskChannelIndex),'Error, you must specify a single MaskChannelIndex for a channel which has valid segmentation!')
end


%% -------- Init ---------- %%
if ~iscell(p.ProcessIndex), p.ProcessIndex={p.ProcessIndex}; end
if ~iscell(p.ChannelIndex), p.ChannelIndex={p.ChannelIndex}; end
if ~iscell(p.OutputName), p.OutputName={p.OutputName}; end

nFrames = movieData.nFrames_;
imSize = movieData.imSize_;

nChan = cellfun(@numel,p.ChannelIndex);
nInput = sum(nChan);
if ~isempty(segProc.maxIndex_)
    maxIndex = segProc.maxIndex_;
else
    maxIndex = 1;
end

%Set up and store the output directories for the window samples.
mkClrDir(p.OutputDirectory)
outFilePaths =cell(numel(p.ProcessIndex),numel(movieData.channels_));
for i=1:numel(p.ProcessIndex)
    iProc = p.ProcessIndex{i};
    if isempty(iProc)
        pString='Raw images - channel ';
    else
        parentOutput = movieData.processes_{iProc}.getDrawableOutput;
        iOutput = strcmp(p.OutputName{i},{parentOutput.var});
        pString=[parentOutput(iOutput).name ' - channel '];
    end
    for j=1:numel(p.ChannelIndex{i})
        iChan = p.ChannelIndex{i}(j);
        outFilePaths{i,iChan} =  [p.OutputDirectory filesep pString num2str(iChan) '.mat'];
    end
end
sampProc.setOutFilePaths(outFilePaths);

%Initialize sample array
sampledFields = {'avg','std','max','min','med','n'};
nFields = numel(sampledFields);
allSamples(nInput,1)=struct();
for j = 1:nFields
    [allSamples.(sampledFields{j})] = deal(nan(maxIndex,1,nFrames));%We add a singleton dimension to be consistent with other sampling processes
end


%% --------- Sampling --------- %%

if ~p.BatchMode && feature('ShowFigureWindows')
    wtBar = waitbar(0,'Please wait, sampling segmented regions...');
else 
    wtBar = -1;
end  

disp('Starting image sampling...');


for iFrame = 1:nFrames
     
    %Load the mask
    currMask = segProc.loadChannelOutput(p.MaskChannelIndex,iFrame);    
        
    
    %Go through each channel and sample it
    stack2sample=zeros([imSize nInput]);

    for i=1:numel(p.ProcessIndex)
        iProc = p.ProcessIndex{i};
        for j=1:numel(p.ChannelIndex{i})
            iChan = p.ChannelIndex{i}(j);
            iInput = sum(nChan(1:i-1))+j;

            if ~isempty(iProc)
                stack2sample(:,:,iInput) = movieData.processes_{iProc}.loadChannelOutput(iChan,iFrame,'output',p.OutputName{i});
            else
                stack2sample(:,:,iInput) = movieData.channels_(iChan).loadImage(iFrame);
            end
        end
    end
    
    %Get the indices present in the current image
    currInd = unique(currMask(currMask(:)>0));
    currSamples = sampleStackSegmentation(currMask,stack2sample,currInd);
    
    assert(numel(currSamples)==nInput);
    
    %Copy these into the whole-movie array    
    for i=1:nInput
        for j = 1:nFields            
            allSamples(i).(sampledFields{j})(currInd,1,iFrame) ...
                = currSamples(i).(sampledFields{j});
        end
    end
       
        
    if ishandle(wtBar) && mod(iFrame,5)
        %Update the waitbar occasionally  to minimize slowdown
        waitbar(iFrame/nFrames,wtBar)
    end
    
end

%% ------- Output ------- %%

if ishandle(wtBar), close(wtBar); end

for i=1:numel(p.ProcessIndex)
    for j=1:numel(p.ChannelIndex{i})
        iChan = p.ChannelIndex{i}(j);
        iInput = sum(nChan(1:i-1))+j;

        samples = allSamples(iInput);  %#ok<NASGU>
        save(outFilePaths{i,iChan},'samples');
    end
end

%Update the movie data, save it
sampProc.setPara(p);%We've stored additional parameters, so add to the process structure.

disp('Finished sampling!')

