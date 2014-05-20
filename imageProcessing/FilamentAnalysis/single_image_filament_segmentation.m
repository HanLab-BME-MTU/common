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
ip.addOptional('pick_channel', 1,@isnum);
ip.parse(filename,Parameter_MD,varargin{:});
pick_channel = ip.Results.pick_channel;
                
%% get the image dir and make a MD just for this image 

filename = GetFullPath(filename);
index_1 = find(filename=='/' | filename=='\');
ROOT_DIR = filename(1:max(index_1));
ROOT_DIR(index_1)=filesep;
file_image_only_name = filename(max(index_1)+1:end);
index_2 = find(file_image_only_name=='.');
file_image_only_name = file_image_only_name(1:max(index_2));

channels_obj_cell = cell(1,numel(Parameter_MD.channels_));

% duplicate channels

for iC = 1 : numel(Parameter_MD.channels_)
    
    % copy the image into a folder for movieData structure
    imagefolder_DIR = [ROOT_DIR,'image',num2str(iC)];
    mkdir(imagefolder_DIR);    
    copyfile(filename,[imagefolder_DIR,filesep,'temp',num2str(iC),'.tif']);
    
    % use this folder for a channel
    channels_obj_cell{iC} = Channel([ROOT_DIR,'image',num2str(iC)]);
    channels_obj_cell{iC}.getImageFileNames();
    channels_obj_cell{iC}.sanityCheck();
    
end

   % build the MD 
    this_MD = MovieData(channels_obj_cell,ROOT_DIR);
    this_MD.sanityCheck();
    this_MD.save();

   
this_MD.addPackage(FilamentAnalysisPackage);

for iPro =  1 : numel(Parameter_MD.processes_)
    if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Thresholding'))
        this_MD.addProcess(ThresholdProcess(this_MD,'funParams',Parameter_MD.processes_{iPro}.funParams_));
        this_MD = thresholdMovie(this_MD,this_MD.processes_{iPro}.funParams_);
        
    else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Mask Refinement'))
            this_MD.addProcess(MaskRefinementProcess(this_MD,'funParams',Parameter_MD.processes_{iPro}.funParams_));
            this_MD = refineMovieMasks(this_MD,this_MD.processes_{iPro}.funParams_);
            
        else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Image Flatten'))
                this_MD.addProcess(ImageFlattenProcess(this_MD,'funParams',Parameter_MD.processes_{iPro}.funParams_));
                this_MD = image_flatten(this_MD);
                
            else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Steerable filtering'))
                    this_MD.addProcess(SteerableFilteringProcess(this_MD,'funParams',Parameter_MD.processes_{iPro}.funParams_));
                    this_MD = steerable_filter_forprocess(this_MD);
                    
                else if(strcmp(Parameter_MD.processes_{iPro}.getName, 'Filament Segmentation'))
                        this_MD.addProcess(FilamentSegmentationProcess(this_MD,'funParams',Parameter_MD.processes_{iPro}.funParams_));
                        this_MD = filament_segmentation(this_MD);
                        
                    end
                end                
            end
        end
    end
end

load([DataOutputDir,'/steerable_vote_1.mat'],...
    'RGB_seg_orient_heat_map', ...
    'MAX_st_res','nms', 'current_seg', ...
    'current_model', 'current_seg_orientation');


save([imagefolder_DIR, filesep, file_image_only_name,'filament_seg_results.mat'],...
    'RGB_seg_orient_heat_map', ...
    'MAX_st_res','nms', 'current_seg', ...
    'current_model', 'current_seg_orientation');
