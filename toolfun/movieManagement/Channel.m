classdef Channel < hgsetget
    %  Class definition of channel class    
    
    properties
        excitationWavelength_       % Excitation wavelength (nm)
        emissionWavelength_         % Emission wavelength (nm)
        exposureTime_               % Exposure time (ms)
        imageType_                  % e.g. Widefield, TIRF, Confocal etc.
        fluorophore_=''               % Fluorophore / Dye (e.g. CFP, Alexa, mCherry etc.)
        
        % ---- Un-used params ---- %
        excitationType_             % Excitation type (e.g. Xenon or Mercury Lamp, Laser, etc)
        neutralDensityFilter_       % Neutral Density Filter
        incidentAngle_              % Incident Angle - for TIRF (degrees)
        filterType_                 % Filter Type
    end
    
    properties(SetAccess=protected)
        psfSigma_                   % Standard deviation of the psf
        channelPath_                % Channel path (directory containing image(s))
        owner_                      % MovieData object which owns this channel
    end
    
    properties(Transient=true)
        displayMethod_  = ImageDisplay; % Method to display the object content
        fileNames_;
    end
    
    methods
        function obj = Channel(channelPath, varargin)
            % Constructor of channel object
            %
            % Input:
            %    channelPath (required) - the absolute path where the channel images are stored
            %
            %    'PropertyName',propertyValue - A string with an valid channel property followed by the
            %    value.
            
            if nargin>0
                obj.channelPath_ = channelPath;
                
                % Construct the Channel object
                nVarargin = numel(varargin);
                if nVarargin > 1 && mod(nVarargin,2)==0
                    for i=1 : 2 : nVarargin-1
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end
        
        %% Set / Get Methods 
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
        
        function set.imageType_(obj, value)
            obj.checkPropertyValue('imageType_',value);
            obj.imageType_=value;
        end
        
        function set.filterType_(obj, value)
            obj.checkPropertyValue('filterType_',value);
            obj.filterType_=value;
        end
        
        function set.fluorophore_(obj, value)
            obj.checkPropertyValue('fluorophore_',value);
            obj.fluorophore_=value;
        end
        
        
        function setFig = edit(obj)
            setFig = channelGUI(obj);
        end
        
        function relocate(obj,oldRootDir,newRootDir)
            % Relocate location of the  channel object
            
            obj.channelPath_=  relocatePath(obj.channelPath_,oldRootDir,newRootDir);
        end
        
        function checkPropertyValue(obj,property, value)
            % Check if a property/value pair can be set up
            
            % Return if unchanged property
            if isequal(obj.(property),value), return; end
            propName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
            
            % Test if the property is writable
            assert(obj.checkProperty(property),'lccb:set:readonly',...
                ['The channel''s ' propName ' has been set previously and cannot be changed!']);
            
            % Test if the supplied value is valid
            assert(obj.checkValue(property,value),'lccb:set:invalid',...
                ['The supplied ' propName ' is invalid!']);
        end
        
                
        function status = checkProperty(obj,property)
            % Returns true/false if the non-empty property is writable
            status = isempty(obj.(property));
            if status, return; end
            
            if strcmp(property,'channelPath_');
                stack = dbstack;
                if any(cellfun(@(x)strcmp(x,'Channel.relocate'),{stack.name})),
                    status  =true;
                end
            end
        end
        
        %---- Sanity Check ----%
        %Verifies that the channel specification is valid, and returns
        %properties of the channel
        
        function [width height nFrames] = sanityCheck(obj,varargin)
            % Check the sanity of the channels
            %
            % Check the validity of each channel and return pixel size and time
            % interval parameters
            
            % Check input
            ip = inputParser;
            ip.addOptional('owner',obj.owner_,@(x) isa(x,'MovieData'));
            ip.parse(varargin{:})
            
            % Set the channel owner
            if isempty(obj.owner_), obj.owner_=ip.Results.owner; end
            assert(isequal(obj.owner_,ip.Results.owner) ||...
                isequal(obj.owner_,ip.Results.owner.parent_),...
                'The channel''s owner is not the movie neither its parent')
            
            if exist(obj.channelPath_, 'file')==2 
                % Using bioformat-tools, get metadata 
                r=bfGetReader(obj.channelPath_,false);
                width = r.getSizeX;
                height = r.getSizeY;
                nFrames = r.getSizeT;
                r.close;
            else
                % Check channel path existence
                assert(logical(exist(obj.channelPath_, 'dir')), ...
                    'Channel path specified is not a valid directory! Please double check the channel path!')
                
                % Check the number of file extensions
                [fileNames nofExt] = imDir(obj.channelPath_,true);
                assert(nofExt~=0,['No proper image files are detected in:'...
                    '\n\n%s\n\nValid image file extension: tif, TIF, STK, bmp, BMP, jpg, JPG.'],obj.channelPath_);
                assert(nofExt==1,['More than one type of image files are found in:'...
                    '\n\n%s\n\nPlease make sure all images are of same type.'],obj.channelPath_);
                nFrames = length(fileNames);
                
                % Check the consistency of image size in current channel
                imInfo = arrayfun(@(x)imfinfo([obj.channelPath_ filesep x.name]),...
                    fileNames, 'UniformOutput', false);
                width = unique(cellfun(@(x)(x.Width), imInfo));
                height = unique(cellfun(@(x)(x.Height), imInfo));
                
                % Check unicity of image sizes
                assert(isscalar(width) && isscalar(height),...
                    ['Image sizes are inconsistent in: \n\n%s\n\n'...
                    'Please make sure all the images have the same size.'],obj.channelPath_);
                obj.fileNames_ = arrayfun(@(x) x.name,fileNames,'UniformOutput',false);
            end
            
            if isempty(obj.psfSigma_) && ~isempty(obj.owner_), obj.calculatePSFSigma(); end
        end
        
        function fileNames = getImageFileNames(obj,iFrame)
            
            if exist(obj.channelPath_, 'file')==2
                % Generate image file names using wavelength if applicable
                [~,channelName]=fileparts(obj.channelPath_);
                if ~isempty(obj.emissionWavelength_)
                    basename = sprintf('%s_w%g_t',channelName ,obj.emissionWavelength_);
                else
                    basename = sprintf('%s_c%d_t',channelName,find(obj.owner_.channels_==obj));
                end
                nFrames=obj.owner_.nFrames_;
                fileNames = arrayfun(@(t) [basename num2str(t, ['%0' num2str(floor(log10(nFrames))+1) '.f']) '.tif'],...
                    1:nFrames,'Unif',false);
            else
                % Channel path is a directory of image files
                if ~isempty(obj.fileNames_)
                    fileNames = obj.fileNames_;
                else
                    fileNames = arrayfun(@(x) x.name,imDir(obj.channelPath_),...
                        'UniformOutput',false);
                end
            end
            if nargin>1, fileNames=fileNames(iFrame); end
        end
        
        function I = loadImage(obj,iFrame)      
            I=zeros([obj.owner_.imSize_ numel(iFrame)]);
            if exist(obj.channelPath_, 'file')==2  
                % Using bioformat tools, get the reader and retrieve dimension order
                loci.common.DebugTools.enableLogging('OFF');
                r=bfGetReader(obj.channelPath_,false);
                metadata=r.getMetadataStore;
                dimensionOrder =char(metadata.getPixelsDimensionOrder(0));
                CZTOrder = dimensionOrder(3:end);
                CZTdimensions = arrayfun(@(x) metadata.(['getPixelsSize' x])(0).getValue,CZTOrder);
                
                % Get channel index, create multi-dimensional plane index and
                % convert into linear plane index
                chanIndex= find(obj.owner_.channels_==obj);
                index(CZTOrder=='C',:)=chanIndex*ones(numel(iFrame),1);
                index(CZTOrder=='Z',:)=1*ones(numel(iFrame),1);
                index(CZTOrder=='T',:)=iFrame;
                iPlane = sub2ind(CZTdimensions,index(1,:),index(2,:),index(3,:));
                
                % Get all requrested planes and close reader
                for i=1:numel(iPlane), I(:,:,i) = bfGetPlane(r,iPlane(i)); end
                r.close;
            else
                % Read images from disk
                fileNames=obj.getImageFileNames(iFrame);
                for i=1:numel(iFrame)
                    I(:,:,i)  = imread([obj.channelPath_ filesep fileNames{i}]);
                end
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
            
            % Input check
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Channel') || numel(x)<=3);
            ip.addRequired('iFrame',@isscalar);
            ip.addParamValue('hAxes',gca,@ishandle);
            ip.KeepUnmatched = true;
            ip.parse(obj,iFrame,varargin{:})
            
            % Initialize output
            if numel(obj)>1, zdim=3; else zdim=1; end
            data = zeros([obj(1).owner_.imSize_ zdim]);
             
            % Fill output
            for iChan=1:numel(obj)
                data(:,:,iChan)=mat2gray(obj(iChan).loadImage(iFrame));
            end
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h = obj(1).displayMethod_.draw(data,'channels','hAxes',ip.Results.hAxes,drawArgs{:});
        end          
    end
    
    methods(Access=protected)
        function calculatePSFSigma(obj)
            % Read parameters for psf sigma calculation
            emissionWavelength=obj.emissionWavelength_*1e-9;
            numAperture=obj.owner_.numAperture_;
            pixelSize=obj.owner_.pixelSize_*1e-9;
            if isempty(emissionWavelength) || isempty(numAperture) || isempty(pixelSize), 
                return; 
            end
            
            obj.psfSigma_ = getGaussianPSFsigma(numAperture,1,pixelSize,emissionWavelength);
%             obj.psfSigma_ =.21*obj.emissionWavelength/(numAperture*pixelSize);

        end
    end
    methods(Static)
        function status=checkValue(property,value)
            % Return true/false if the value for a given property is valid
            
            % Parse input
            ip = inputParser;
            ip.addRequired('property',@(x) ischar(x) || iscell(x));
            ip.parse(property);
            if iscell(property)
                ip.addRequired('value',@(x) iscell(x)&&isequal(size(x),size(property)));
                ip.parse(property,value);
                status=cellfun(@(x,y) Channel.checkValue(x,y),property,value);
                return
            end
            
            % Get validator for single property
            validator=Channel.getPropertyValidator(property);
            propName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
            assert(~isempty(validator),['No validator defined for property ' propName]);
            
            % Return result of validation
            status = isempty(value) || validator(value);
        end
        
        function validator = getPropertyValidator(property)
            switch property
                case {'emissionWavelength_','excitationWavelength_'}
                    validator=@(x) isscalar(x) && x>=300 && x<=800;
                case 'exposureTime_'
                    validator=@(x) isscalar(x) && x>0;
                case {'excitationType_','notes_','channelPath_','filterType_'}
                    validator=@ischar;
                case 'imageType_'
                    validator = @(x) ischar(x) && ismember(x,Channel.getImagingModes);
                case {'fluorophore_'}
                    validator= @(x) ischar(x) && ismember(x,Channel.getFluorophores);
                otherwise
                    validator=[];
            end
            
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
