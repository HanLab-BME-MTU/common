classdef Channel < hgsetget
    %  Class definition of channel class
    properties (SetAccess = protected)
        psfSigma_                  % standard deviation of the psf
    end
        
    properties 
        
        % ---- Used Image Parameters ---- %
                excitationWavelength_       % Excitation wavelength (nm)
        emissionWavelength_         % Emission wavelength (nm)
        exposureTime_               % Exposure time (ms)
        imageType_                  % e.g. Widefield, TIRF, Confocal etc.
        
        % ---- Un-used params ---- %
        
        excitationType_             % Excitation type (e.g. Xenon or Mercury Lamp, Laser, etc)
        neutralDensityFilter_       % Neutral Density Filter
        incidentAngle_              % Incident Angle - for TIRF (degrees)
        filterType_                 % Filter Type
        fluorophore_=''               % Fluorophore / Dye (e.g. CFP, Alexa, mCherry etc.)  
        
    end
    
    properties(SetAccess=protected) 
        % ---- Object Params ---- %
        channelPath_                % Channel path (directory containing image(s))
        owner_                      % MovieData object which owns this channel 
    end
    
    properties(Transient=true)
        displayMethod_  = ImageDisplay;
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

        function relocate(obj,oldRootDir,newRootDir)

            % Relocate channel path
            obj.channelPath_=  relocatePath(obj.channelPath_,oldRootDir,newRootDir);
            
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
                if strcmp(stack(3).name,'Channel.relocate'), propertyCheck=1; end
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
        
        function [width height nFrames] = sanityCheck(obj,varargin)
            % Check the validity of each channel and return pixel size and time
            % interval parameters
            
            % Check input
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Channel'));
            ip.addOptional('owner',obj.owner_,@(x) isa(x,'MovieData'));
            ip.parse(obj,varargin{:})
            obj.owner_=ip.Results.owner;
            
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
            imInfo = arrayfun(@(x)imfinfo([obj.channelPath_ filesep x.name]),...
                fileNames, 'UniformOutput', false);
            width = unique(cellfun(@(x)(x.Width), imInfo));
            height = unique(cellfun(@(x)(x.Height), imInfo));
            
            % Exception: Image sizes are inconsistent in the
            % current channel.
            assert(isscalar(width) && isscalar(height),...
                ['Image sizes are inconsistent in: \n\n%s\n\n'...
                'Please make sure all the images have the same size.'],obj.channelPath_);
            
            if isempty(obj.psfSigma_) && ~isempty(obj.owner_)
                obj.calculatePSFSigma(obj.owner_.numAperture_,obj.owner_.pixelSize_);
            end
               
        end
        
        function fileNames = getImageFileNames(obj,iFrame)
            fileNames = arrayfun(@(x) x.name,imDir(obj.channelPath_),...
                'UniformOutput',false);
            if nargin>1
                fileNames=fileNames{iFrame};
            end
        end
        
        function color = getColor(obj)

            if ~isempty(obj.emissionWavelength_),
                color = wavelength2rgb(obj.emissionWavelength_*1e-9);
            else
                color =[1 1 1]; % Set to grayscale by default
            end
        end
        
        function h = draw(obj,iFrame,varargin)
           
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Channel'));
            ip.addRequired('iFrame',@isscalar);
            if numel(obj)==1
                ip.addParamValue('color',obj.getColor,@(x) isequal(size(x),[1 3]));
            end
            ip.addParamValue('hAxes',gca,@ishandle);
            ip.parse(obj,iFrame,varargin{:})
            
            if numel(obj)>3, error('Max. 3 channels in RGB mode.'); end
            if numel(obj)>1
                data = zeros([obj(1).owner_.imSize_ 3]);
                for iChan=1:numel(obj)
                    imageName = obj(iChan).getImageFileNames(iFrame);
                    rawData = double(imread([obj(iChan).channelPath_ filesep imageName]));
                    data(:,:,iChan)=rawData/max(rawData(:));
                end
            else
                imageName = obj.getImageFileNames(iFrame);
                rawData = double(imread([obj.channelPath_ filesep imageName]));
                data=repmat(rawData/max(rawData(:)),[1 1 3]);
                color=ip.Results.color;
                for i=1:3
                    data(:,:,i)=data(:,:,i)*color(i);
                end
            end
%             
%             frame = zeros(obj.Im,nx,3);
%             idxRGB = getRGBindex(data.markers);
%             for c = 1:nCh
%                 frame(:,:,idxRGB(c)) = scaleContrast(double(imread(data.framePaths{c}{frameIdx})), ip.Results.iRange{c});
%             end
%             imageName = obj.getImageFileNames(iFrame);
%             image = double(imread([obj.channelPath_ filesep imageName]));
            
            h = obj(1).displayMethod_.draw(data,'channels','hAxes',ip.Results.hAxes);
        end
        
        
    end
    
    methods(Access=protected)
        function calculatePSFSigma(obj,numAperture,pixelSize)
            if isempty(obj.emissionWavelength_); return; end
            if isempty(numAperture) || isempty(pixelSize), return; end
            if strcmp(obj.imageType_,'Widefield')
                obj.psfSigma_ =.21*obj.emissionWavelength_/(numAperture*pixelSize);
            else
                obj.psfSigma_ = getGaussianPSFsigma(numAperture,1,...
                    pixelSize*1e-9,obj.emissionWavelength_*1e-9);
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
            
            if iscell(property)
                checkValue=cellfun(@(x,y) Channel.checkValue(x,y),property,value);
                return
            end
            
            switch property                
                case {'emissionWavelength_','excitationWavelength_'}
                    checkTest=@(x) isnumeric(x) && x>=300 && x<=800;
                case 'exposureTime_'
                    checkTest=@(x) isnumeric(x) && x>0;
                case {'excitationType_','notes_','channelPath_','filterType_'}
                    checkTest=@(x) ischar(x);
                case 'imageType_'
                    checkTest = @(x) any(strcmpi(x,Channel.getImagingModes));
                case {'fluorophore_'}
                    checkTest= @(x)  any(strcmpi(x,Channel.getFluorophores));
                case {'owner_'}
                    checkTest= @(x) isa(x,'MovieData');
            end
            checkValue = isempty(value) || checkTest(value);
        end
        
        function modes=getImagingModes()
            modes={'Widefield';'TIRF';'Confocal'};
        end
        
        function fluorophores=getFluorophores()
            fluorPropStruct= getFluorPropStruct();
            fluorophores={fluorPropStruct.name};
        end
        
    end
end
