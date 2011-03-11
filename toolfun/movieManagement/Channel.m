classdef Channel < handle
%  Class definition of channel class    

    properties(SetAccess = protected, GetAccess = public)
        
        % ---- Used Image Parameters ---- %
        
        channelPath_                % Channel path (directory containing image(s))
        excitationWavelength_       % Excitation wavelength (nm)
        emissionWavelength_         % Emission wavelength (nm)
        exposureTime_               % Exposure time (ms)
        imageType_                  % e.g. Widefield, TIRF, Confocal etc.
        
        % ---- Un-used params ---- %
        
        excitationType_             % Excitation type (e.g. Xenon or Mercury Lamp, Laser, etc)
        neutralDensityFilter_       % Neutral Density Filter
        incidentAngle_              % Incident Angle - for TIRF (degrees)
        filterType_                 % Filter Type
        fluorophore_                % Fluorophore / Dye (e.g. CFP, Alexa, mCherry etc.)
        
        % ---- Object Params ---- %
        owner_                      %MovieData object which owns this channel
        
        
    end
    
    methods (Access = public)
        
        function obj = Channel (channelPath, excitationWavelength, ...
                            emissionWavelength, exposureTime, excitationType, ...
                            neutralDensityFilter, incidentAngle, filterType, fluorophore)
            if nargin > 0
                obj.channelPath_ = channelPath;
            end

            if nargin > 1
                obj.excitationWavelength_ = excitationWavelength;
            end
            if nargin > 2
                obj.emissionWavelength_ = emissionWavelength;
            end
            if nargin > 3
                obj.exposureTime_ = exposureTime;
            end
            if nargin > 4
                obj.excitationType_ = excitationType;
            end
            if nargin > 5
                obj.neutralDensityFilter_ = neutralDensityFilter;
            end
            if nargin > 6
                obj.incidentAngle_ = incidentAngle;
            end         
            if nargin > 7
                obj.filterType_ = filterType;
            end   
            if nargin > 8
                obj.fluorophore_ = fluorophore;
            end                   

        end
        
        function setChannelPath(obj,chanPath)
            if isempty(obj.channelPath_)
                obj.channelPath_ = chanPath;
            else
                % Check the calling function
                stack = dbstack;
                if  strcmp(stack(2).name,'MovieData.relocateMovieData')
                    obj.channelPath_ = chanPath;
                else
                    error('This channel''s path has already been set and cannot be changed!');
                end
            end
        end
        function setExcitationWavelength(obj, data)
            obj.excitationWavelength_ = data;
        end
        
        function setEmissionWavelength(obj, data)
            obj.emissionWavelength_ = data;
        end
        
        function setExposureTime(obj, data)
            obj.exposureTime_ = data;
        end  

        function setExcitationType(obj, data)
            obj.excitationType_ = data;
        end  

        function setNeutralDensityFilter(obj, data)
            obj.neutralDensityFilter_ = data;
        end  
        
        function setIncidentAngle(obj, data)
            obj.incidentAngle_ = data;
        end  
        
        function setFilterType(obj, data)
            obj.filterType_ = data;
        end  
        
        function setFluorophore(obj, data)
            obj.fluorophore_ = data;
        end          
        
        function setOwner(obj,owner)
            if isempty(obj.owner_)
                if isa(owner,'MovieData')
                    obj.owner_ = owner;
                else
                    error('The Channel can only be owned by a MovieData object!')
                end
                
            else
                error('This channel already has an owner, and this property cannot be changed!');
            end
        end
        function [width height nFrames nSlices] = sanityCheck(obj)
        % Check the validity of each channel and return pixel size and time
        % interval parameters
            
            % Exception: channel path does not exist
            assert(logical(exist(obj.channelPath_, 'dir')), ...
                    'One or more channel paths do not exist. Please double check the channel paths.')
            
            % Check the number of file extensions
            [fileNames nofExt] = imDir(obj.channelPath_,true);
            switch nofExt
                case 0
                    % Exception: No proper image files are detected
                    error('No proper image files are detected in:\n\n%s\n\nValid image file extension: tif, TIF, STK, bmp, BMP, jpg, JPG.',obj.channelPath_);
                                                  
                case 1
                    nFrames = length(fileNames); 
                    
                otherwise
                    % Exception: More than one type of image 
                    % files are in the current specific channel
                    error('More than one type of image files are found in:\n\n%s\n\nPlease make sure all images are of same type.', obj.channelPath_); 
            end
            
            % Check the consistency of image size in current channel            
            imInfo = arrayfun(@(x)imfinfo([obj.channelPath_ filesep x.name]), fileNames, 'UniformOutput', false);
            imSize2(1,:) = cellfun(@(x)(x.Width), imInfo, 'UniformOutput', true);
            imSize2(2,:) = cellfun(@(x)(x.Height), imInfo, 'UniformOutput', true);
            
            % Exception: Image sizes are inconsistent in the
            % current channel. 
            
            assert(max(imSize2(1,:))==min(imSize2(1,:)) && ...
                    max(imSize2(2,:))==min(imSize2(2,:)), ...
                    'Image sizes are inconsistent in: \n\n%s\n\nPlease make sure all the images have the same size.',obj.channelPath_)

            width = imSize2(1);
            height = imSize2(2);            
            
            %If it's a 3D movie, make sure the number of slices is correct
            if ~isempty(obj.owner_) && isa(obj.owner_,'MovieData3D')

                imNames = imDir(obj.channelPath_);
                
                %I guess we have to read an image to check this...?
                tmp = size(stackRead([obj.channelPath_ filesep imNames(1).name])); %Just check the first image... Lazy but tooooo slow to check all!! Need to figure out how to check without loading....
                if numel(tmp) ~= 3 || min(tmp) < 2
                    error('The images in the specified image folder are not 3D stacks!')
                end
                nSlices = tmp(3);
                if nSlices ~= obj.owner_.nSlices_
                    error('The number of z sections found does not match that specified by owner MovieData!')
                end
                
            else
                nSlices = [];
            end
        end
        
    end
    
end