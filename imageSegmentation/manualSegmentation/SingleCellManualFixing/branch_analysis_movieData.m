function BA_output_cell = branch_analysis_movieData(MD,half_size,min_branch_size_Threshold,figure_flag)
% function to do branch analysis for a whole movieData
% Liya Ding, March, 2014
%
% Input:
%   ML:     The movieList object loaded before running this function


if(nargin<2)
    half_size=150;
end

if(nargin<3)
    min_branch_size_Threshold=100;
end

if(nargin<4)
    figure_flag=0;
end

BA_output_cell = cell(1,1);

% the number of channels
nChannel = length(MD.channels_);

ROOT_DIR = MD.outputDirectory_;

for iChannel = 1 : nChannel
    for iCell = 1 : 20        
        display(['iChannel:', num2str(iChannel),', iCell:', num2str(iCell)]);
        
        FilamentAnalysisPackage_complete_frames_file_name = [ROOT_DIR,'\','FilamentAnalysisPackage','\completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),'\completedFrames.mat'];
        SegmentationPackage_complete_frames_file_name = [ROOT_DIR,'\','SegmentationPackage','\completedFramesChannel',num2str(iChannel),'Cell',num2str(iCell),'\completedFrames.mat'];
        PackageName=[];
        
        if(exist(FilamentAnalysisPackage_complete_frames_file_name,'file'))
            PackageName = 'FilamentAnalysisPackage';
        end
        
        if(exist(SegmentationPackage_complete_frames_file_name,'file'))
            PackageName = 'SegmentationPackage';
        end
        
        if(isempty(PackageName))
            continue;
        end
        
        % the folder name if there is marking
        truthPath = [ROOT_DIR,'\',PackageName,'\FixedChannel',num2str(iChannel),'Cell',num2str(iCell)];
        % check if this folder exist
        if(exist(truthPath,'dir'))
            % if it exist, try to do the branch analysis
            BA_output = branch_analysis_marked_cell(MD, iChannel, iCell,half_size,figure_flag);
            % for the output, ignore the channel number since
            % there is only one channel marked.
            if(~isempty(BA_output))
                BA_output_cell{iChannel, iCell} = BA_output;            
            end
        end
    end
end


save([ROOT_DIR,'\movieData_BA_output.mat'],'BA_output_cell');

