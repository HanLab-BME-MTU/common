function MLorMDOut=addAnalysisFolder(MLorMD,currentAnalysisRoot,newAnalysisRoot,varargin)
% This function copies a movieData object (or ML) and creates a new analysis folder.
% It never re-write the orignal MD or ML file. It refers to the same
% channels.
%
% Optionnaly the channel can be relocated to using the options
% oldRawDataRoot and newRawDataRoot.
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched=true;
ip.addRequired('MLorMD');
ip.addRequired('currentAnalysisRoot',@ischar);
ip.addRequired('newAnalysisRoot',@ischar);
ip.addOptional('oldRawDataRoot','',@ischar);
ip.addOptional('newRawDataRoot','',@ischar);
ip.addParamValue('recursive',true,@islogical);
ip.parse(MLorMD,currentAnalysisRoot,newAnalysisRoot,varargin{:});
oldRawDataRoot=ip.Results.oldRawDataRoot;
newRawDataRoot=ip.Results.newRawDataRoot;

if( isa(MLorMD,'MovieData'))
    MLorMDOut=addAnalysisFolderMD(MLorMD,currentAnalysisRoot,newAnalysisRoot,varargin{:});
else    
    ML=MLorMD;
    MDs=cell(1,length(ML.movieDataFile_));
    for i=1:length(ML.movieDataFile_)
        MD=[];
        if(~isempty(oldRawDataRoot))
            MD=MovieData.loadMatFile(relocatePath(ML.movieDataFile_{i},oldRawDataRoot,newRawDataRoot));
        else
            MD=MovieData.loadMatFile(ML.movieDataFile_{i});
        end
        MDs{i}=addAnalysisFolderMD(MD,currentAnalysisRoot,newAnalysisRoot,varargin{:});
        
    end
    mkdir([newAnalysisRoot filesep 'analysis']);
    MLorMDOut=MovieList(MDs,[newAnalysisRoot filesep 'analysis'],'movieListFileName_',ML.movieListFileName_,'movieListPath_',[newAnalysisRoot filesep 'analysis']);
    MLorMDOut.save();
end

function MD=addAnalysisFolderMD(MD,currentAnalysisRoot,newAnalysisRoot,varargin)
    ip = inputParser;
    ip.CaseSensitive = false;
    ip.KeepUnmatched=true;
    ip.addRequired('MD');
    ip.addRequired('currentAnalysisRoot',@ischar);
    ip.addRequired('newAnalysisRoot',@ischar);
    ip.addOptional('oldRawDataRoot','',@ischar);
    ip.addOptional('newRawDataRoot','',@ischar);
    ip.parse(MD,currentAnalysisRoot,newAnalysisRoot,varargin{:});
    
    oldRawDataRoot=ip.Results.oldRawDataRoot;
    newRawDataRoot=ip.Results.newRawDataRoot;
    
    MDAnalysisPath=relocatePath(MD.outputDirectory_, currentAnalysisRoot,  newAnalysisRoot);
    mkdir(MDAnalysisPath);
    %MD.sanityCheck(MDAnalysisPath,[sprintf(namePattern,i) '.mat'],false);
    % , filesep sprintf(namePattern,i)
    
    for c=1:length(MD.channels_);
        MD.getChannel(c).relocate(oldRawDataRoot,newRawDataRoot);
    end
    mkdir(MDAnalysisPath);
    MD.relocate(MD.outputDirectory_,MDAnalysisPath,false);
    MD.save();
