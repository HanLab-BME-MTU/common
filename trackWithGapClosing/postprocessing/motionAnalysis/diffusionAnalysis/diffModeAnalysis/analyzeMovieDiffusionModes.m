function analyzeMovieDiffusionModes(movieDataOrProcess,varargin)
% analyzeMovieDiffusionModes performs diffusion mode analysis on single particle tracks
%
% SYNOPSIS analyzeMovieDiffusionModes(movieDataOrProcess,varargin)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below:
%
%   Possible Parameter Structure Field Names:
%       ('FieldName' -> possible values)
%
%       ('OutputDirectory' -> Character string)
%       A character string specifying the directory where to save 
%       the motion analysis results.
%
%       ('ChannelIndex' -> Positive integer scalar or vector)
%       The integer index of the channel(s) containing tracks to be analyzed.
%
%       PARAMETERS FOR DIFFUSION MODE DECOMPOSITION:   
%
%       ('probDim' -> Positive integer) 
%       Problem dimensionality. Default is 2.
% 
%       ('minLength' -> Positive integer) 
%       Minimum length of a track to be included in analysis. 
%       Default: 5 frames.
% 
%       ('alpha' -> Positive double) 
%       Alpha-value for the statistical test to determine number of modes. 
%       Default: 0.01.
%
%       ('showPlot' -> Boolean)
%       1 to plot the histogram of square displacements and fitted exponentials.
%       0 to not plot anything. Default: 1.
%
%       ('maxNumMode' -> Positive integer)
%       Upper limit on the number of modes. Default: 10.            
%
%       ('binStrategy' -> Positive integer)
%       Binning strategy for calculating the cumulative histogram. 1 for 
%       using "histogram" and 2 for using the data directly. Default: 2.
%
%       ('plotName' -> Character string)
%       The title of the plotted figure. Default: 'Figure'.
%
%       ('subSampSize' -> Positive integer)
%       Size of subsample to use in mode decomposition. In this case, the
%       original data are subsampled many times, each time with the
%       specified subSampSize. The output is then the diffusion mode
%       decomposition result for each of the subsamples. Enter [] if no
%       sub-sampling. Default: [].
%
%       ('doControl' -> Boolean)
%       1 in order to do mono-exponential control, 0 otherwise. 
%       Default: 1.
%
%       ('forceDecompose' -> Boolean)
%       1 in order to force run diffusion mode decomposition, 0 to
%       use old decomposition results if available.  
%       Default: 1.
%
%       PARAMETERS FOR TRACK DIFFUSION MODE ASSIGNMENT
%
%       ('diffModeDividerStruct' -> Structure)
%       See function trackDiffModeAnalysis for structure details.       
%       Default: [], in which case the tracks are not classified into modes
%       but only their diffusion coefficient is output.
%
% Khuloud Jaqaman, August 2017

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @isProcessOrMovieData);
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess,varargin{:});
paramsIn=ip.Results.paramsIn;

% Get the DiffusionModeAnalysisProcess and create it if it does not exist
[movieData, postProc] = getOwnerAndProcess(movieDataOrProcess,'DiffusionModeAnalysisProcess',true);

%Parse input, store in parameter structure
p = parseProcessParams(postProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check tracking process first
iTrackProc = movieData.getProcessIndex('TrackingProcess',1,1);

assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to post-processing!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing tracking output! Please apply tracking before ' ...
    'running  post-processing!']);

% Check whether there is a cell mask to use
iProcMask = movieData.getProcessIndex('ImportCellMaskProcess',1,0); %cell mask
if ~isempty(iProcMask)
    mask = imread(fullfile(movieData.processes_{iProcMask}.funParams_.OutputDirectory,'cellMask_channel_1.tif'));
else
    mask = [];
end

% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = trackProc.outFilePaths_{1,i};
end
postProc.setInFilePaths(inFilePaths);

%variables related to diffusion mode decomposition
newDecompose = ones(numel(movieData.channels_),1);
diffModeDecomposeOld = cell(numel(movieData.channels_),1);

% Set up the output file
%store old diffusion mode decomposition if relevant
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
    if ~p.forceDecompose && exist(outFilePaths{1,i})~=0
        tmp = load(outFilePaths{1,i});
        if isfield(tmp,'diffModeDecomposition')
            diffModeDecomposeOld{i} = tmp.diffModeDecomposition;
        else
            diffModeDecomposeOld{i} = tmp.diffModeAnalysisRes;
        end
        newDecompose(i) = 0;
    end
end
mkClrDir(p.OutputDirectory);
postProc.setOutFilePaths(outFilePaths);

%% --------------- Diffusion Mode Analysis ---------------%%% 

disp('Starting diffusion mode analysis ...')

for i = p.ChannelIndex
    
    tracks = trackProc.loadChannelOutput(i);
    
    %     try
    %         if postProc.funParams_.driftCorrect ==1
    %             [tracks] = scriptCorrectImageDrift(tracks,movieData);
    %         end
    %     catch
    %
    %     end
    
    %diffusion mode decomposition from frame-to-frame displacements
    if newDecompose(i)==1
        [modeParam,numMode,modeParamControl,numModeControl] = ...
            getDiffModes(tracks,p.minLength,p.alpha,p.showPlot,p.maxNumMode,...
            p.binStrategy,p.plotName,p.subSampSize,p.doControl,mask);
        diffModeDecomposition.modeParam = modeParam;
        diffModeDecomposition.numMode = numMode;
        diffModeDecomposition.modeParamControl = modeParamControl;
        diffModeDecomposition.numModeControl = numModeControl;
    else
        diffModeDecomposition = diffModeDecomposeOld{i};
    end
    
    %diffusion mode assignment per track
    diffModeAnalysisRes = trackDiffModeAnalysis(tracks,p.diffModeDividerStruct);
    
    tracks = trackProc.loadChannelOutput(i);
    for j = 1:numel(tracks)
        tracks(j).classification = diffModeAnalysisRes(j).diffMode;
    end
    
    % save each projData in its own directory
    save(outFilePaths{1,i},'diffModeAnalysisRes', 'diffModeDecomposition','tracks')
    
end

disp('Finished diffusion mode analysis!')

