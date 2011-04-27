classdef Channel < hgsetget
    %  Class definition of channel class
    
    properties 
        
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
    end
    
    properties(SetAccess=protected) 
        % ---- Object Params ---- %
        
        owner_                      % MovieData object which owns this channel 
    end
    
    methods
                
        function obj = Channel(channelPath, varargin)
            % Constructor of channel object
            %
            % INPUT  
            %
            %    channelPath (required) - the absolute path where the channel images are stored
            %
            %    'PropertyName',propertyValue - A string with an option name followed by the
            %    value for that option.
            %    Possible Option Names are the Channel fieldnames
            %
            %
            % OUTPUT
            %
            %    obj - an object of class Channel
            %

            obj.channelPath_ = channelPath;

            % Construct the Channel object
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
            
        end
        
        % ------- Set / Get Methods ----- %
        
        function set.channelPath_(obj,value)
            obj.checkPropertyValue('channelPath_',value);
            obj.channelPath_=value;
        end

        function set.excitationWavelength_(obj, value)
            obj.checkPropertyValue('excitationWavelength_',value);
            obj.excitationWavelength_=value;
        end
        
        function set.emissionWavelength_(obj, value)
            obj.checkPropertyValue('emissionWavelength_',value);
            obj.emissionWavelength_=value;
        end
        
        function set.exposureTime_(obj, value)
            obj.checkPropertyValue('exposureTime_',value);
            obj.exposureTime_=value;
        end
        
        function set.excitationType_(obj, value)
            obj.checkPropertyValue('excitationType_',value);
            obj.excitationType_=value;
        end
        
        function set.neutralDensityFilter_(obj, value)
            obj.checkPropertyValue('neutralDensityFilter_',value);
            obj.neutralDensityFilter_=value;
        end
        
        function set.incidentAngle_(obj, value)
            obj.checkPropertyValue('incidentAngle_',value);
            obj.incidentAngle_=value;
        end
        
        function set.filterType_(obj, value)
            obj.checkPropertyValue('filterType_',value);
            obj.filterType_=value;
        end
        
        function set.fluorophore_(obj, value)
            obj.checkPropertyValue('fluorophore_',value);
            obj.fluorophore_=value;
        end
        
        function set.owner_(obj,value)
            obj.checkPropertyValue('owner_',value);
            obj.owner_=value;
        end

        function setFig = edit(obj)
            setFig = channelGUI(obj);
        end

        function checkPropertyValue(obj,property, value)
            % Check if a property/value pair can be set up
            % 
            % Returns an error if either the property is unchangeable or
            % the value is invalid.
            %
            % INPUT:
            %    property - a valid Channel property name (string)
            %
            %    value - the property value to be checked
            %
            
            % Test if the property is writable
            propertyCheck =0;
            if strcmp(property,{'notes_'}), propertyCheck=1;
            elseif isempty(obj.(property)), propertyCheck=1; 
            elseif isequal(obj.(property),value), return;               
            elseif strcmp(property,'channelPath_')
                % Allow relocation of channelPath_
                stack = dbstack;
                if strcmp(stack(3).name,'MovieData.relocate'), propertyCheck=1; end
            end
            
            if ~propertyCheck
                propertyName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
                error(['This channel''s ' propertyName ' has been set previously and cannot be changed!']);
            end
            
            % Test if the value is valid
            valueCheck=obj.checkValue(property,value);
            if ~valueCheck
                propertyName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
                error(['The supplied ' propertyName ' is invalid!']);
            end
        end
        
        %---- Sanity Check ----%
        %Verifies that the channel specification is valid, and returns
        %properties of the channel
        
        function [width height nFrames] = sanityCheck(obj,owner)
            % Check the validity of each channel and return pixel size and time
            % interval parameters
            
            % Exception: channel path does not exist
            assert(logical(exist(obj.channelPath_, 'dir')), ...
                'Channel path specified is not a valid directory! Please double check the channel path!')
            
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
            
            if nargin>1
                if isempty(obj.owner_), 
                    obj.owner_=owner; 
                else
                    assert(obj.owner_==owner,'Channel object can only be owned by one MovieData object');
                end
                
            end
                
        end
        
    end
    
    methods(Static)
        function checkValue=checkValue(property,value)
            % Test the validity of a property value
            %
            % Declared as a static method so that the method can be called
            % without instantiating a Channel object. Should be called
            % by any generic set method.
            %
            % INPUT:
            %    property - a Channel property name (string)
            %
            %    value - the property value to be checked
            %
            % OUTPUT:
            %    checkValue - a boolean containing the result of the test
            
            switch property                
                case {'emissionWavelength_','excitationWavelength_'}
                    checkTest=@(x) isnumeric(x) && x>=300 && x<=800;
                case 'exposureTime_'
                    checkTest=@(x) isnumeric(x) && x>0;
                case {'imageType_','excitationType_','fluorophore_','notes_','channelPath_'}
                    checkTest=@(x) ischar(x);
                case {'owner_'}
                    checkTest= @(x) isa(x,'MovieData');
            end
            checkValue = isempty(value) || checkTest(value);
        end
    end
end
