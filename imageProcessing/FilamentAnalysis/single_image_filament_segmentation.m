function [current_seg,current_seg_orientation, current_model,RGB_seg_orient_heat_map,MAX_st_res,nms] ...
    = single_image_filament_segmentation(filename, Parameter_MD, varargin)
% Function of single image filament segmentation with input MD from other
%               successfully segmented movie for the parameters in MD

% Input:      filename:        the filename of the image to be segmented
%             Parameter_MD:    a loaded MD with good segmentation parameters
%             pick_channel:    optional input, which channel to use in case there are more
%                               than one channel in the MD, default 1
% Output:     current_model: the filament-by-filament model
%             current_seg: the black-white segmentation results
%             current_seg_orientation: the orienation of segmented pixels
%             RGB_seg_orient_heat_map: heat map for display
%             MAX_st_res: is the steerable filtering result
%             nms: the non-maximum-surpress ST

% Created 05 2014 by Liya Ding, Matlab R2012b

ip = inputParser;
ip.addRequired('filename',@ischar);
ip.addRequired('Parameter_MD',@(x) isa(x,'MovieData'));
ip.addOptional('pick_channel', 1,@isnumeric);
ip.addOptional('keep_steps', 0,@islogic);
ip.parse(filename,Parameter_MD,varargin{:});
pick_channel = ip.Results.pick_channel;
keep_steps = ip.Results.keep_steps;

%% get the image dir and make a MD just for this image

filename = GetFullPath(filename);
index_1 = find(filename=='/' | filename=='\');
ROOT_DIR = filename(1:max(index_1));
ROOT_DIR(index_1)=filesep;

file_image_full_name = filename(max(index_1)+1:end);
index_2 = find(file_image_full_name=='.');
file_image_only_name = file_image_full_name(1:max(index_2)-1);

channels_obj_cell = [];
channels_obj_all=[];
% duplicate channels

for iC = 1 : numel(Parameter_MD.channels_)
    
    % copy the image into a folder for movieData structure
    imagefolder_DIR = [ROOT_DIR,'single_run_image',num2str(iC)];
    if (~exist(imagefolder_DIR,'dir'))
        mkdir(imagefolder_DIR);
    end
    copyfile(filename,[imagefolder_DIR,filesep,'temp',num2str(iC),'.tif']);
    
    % use this folder for a channel
    channels_obj_cell = Channel([ROOT_DIR,'single_run_image',num2str(iC)]);
    channels_obj_cell.getImageFileNames();
    channels_obj_cell.sanityCheck();
    
    
    channels_obj_all = [channels_obj_all channels_obj_cell];
end

% build the MD
this_MD = MovieData(channels_obj_all,ROOT_DIR);
this_MD.setPath(ROOT_DIR);
this_MD.setFilename('movieData.mat')
this_MD.sanityCheck();


this_MD.addPackage(FilamentAnalysisPackage);
to_delete_folders = cell(1,1);
folder_count=0;

for iPro =  1 : numel(Parameter_MD.processes_)
    if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Thresholding'))
        given_Params = Parameter_MD.processes_{iPro}.funParams_;
        given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'thres'];
        
        given_Params.ChannelIndex = pick_channel;
        this_MD.addProcess(ThresholdProcess(this_MD,'funParams',given_Params));
        this_MD = thresholdMovie(this_MD,this_MD.processes_{iPro}.funParams_);
        
        folder_count = folder_count+1;
        to_delete_folders{folder_count}= given_Params.OutputDirectory ;
        
    else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Mask Refinement'))
            given_Params = Parameter_MD.processes_{iPro}.funParams_;
            given_Params.OutputDirectory = [this_MD.outputDirectory_,filesep,'maskrefine'];
            given_Params.ChannelIndex = pick_channel;
            this_MD.addProcess(MaskRefinementProcess(this_MD,'funParams',given_Params));
            this_MD = refineMovieMasks(this_MD,this_MD.processes_{iPro}.funParams_);
            
            folder_count = folder_count+1;
            to_delete_folders{folder_count}= given_Params.OutputDirectory ;
            
        else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Image Flatten'))
                given_Params = Parameter_MD.processes_{iPro}.funParams_;
                given_Params.ChannelIndex = pick_channel;
                this_MD.addProcess(ImageFlattenProcess(this_MD));
                this_MD = image_flatten(this_MD,given_Params);
                
                folder_count = folder_count+1;
                to_delete_folders{folder_count}= [this_MD.outputDirectory_,filesep,'ImageFlatten'];
        
                
            else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Steerable filtering'))
                    given_Params = Parameter_MD.processes_{iPro}.funParams_;
                    given_Params.ChannelIndex = pick_channel;
                    this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',given_Params));
                    this_MD = steerable_filter_forprocess(this_MD,given_Params);
                    
                    STOutputDir = [this_MD.processes_{iPro}.outFilePaths_{pick_channel}];
                    
                    load([STOutputDir,'/steerable_temp',num2str(pick_channel),'.tif.mat'],...
                        'MAX_st_res','nms');
                    
                    folder_count = folder_count+1;
                    to_delete_folders{folder_count}= [this_MD.outputDirectory_, filesep, 'SteerableFiltering'];
                        
                else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                        given_Params = Parameter_MD.processes_{iPro}.funParams_;
                        given_Params.ChannelIndex = pick_channel;
                        this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',given_Params));
                        this_MD = filament_segmentation(this_MD,given_Params);
                        
                        
                        DataOutputDir = [this_MD.processes_{iPro}.outFilePaths_{pick_channel},'/DataOutput'];
                        
                        load([DataOutputDir,'/steerable_vote_temp',num2str(pick_channel),'.tif.mat'],...
                            'RGB_seg_orient_heat_map', ...
                            'current_seg', ...
                            'current_model', 'current_seg_orientation');
                        
                        folder_count = folder_count+1;
                        to_delete_folders{folder_count}= [this_MD.outputDirectory_, filesep, 'FilamentSegmentation'];
             
                        
                    end
                end
            end
        end
    end
end


save([ROOT_DIR, filesep, file_image_only_name,'_filament_seg_results.mat'],...
    'RGB_seg_orient_heat_map', ...
    'MAX_st_res','nms', 'current_seg', ...
    'current_model', 'current_seg_orientation');

if(keep_steps==0)
    
    for iF =1 : folder_count
    rmdir(to_delete_folders{iF},'s');
    end
end

