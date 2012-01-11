classdef  MovieData < MovieObject
    % Concrete implementation of MovieObject for a single movie
    
    properties (SetAccess = protected)
        channels_ = [];         % Channel object array
        nFrames_                % Number of frames
        imSize_                 % Image size 1x2 array[height width]
    end
    
    properties
        movieDataPath_          % The path where the movie data is saved
        movieDataFileName_      % The name under which the movie data is saved
        pixelSize_              % Pipxel size (nm)
        timeInterval_           % Time interval (s)
        numAperture_            % Numerical Aperture
        camBitdepth_            % Camera Bit-depth
        eventTimes_             % Time of movie events
        
        % ---- Un-used params ----
        
        magnification_
        binning_
        
    end
    
    methods
        %% Constructor
        function obj = MovieData(channels,outputDirectory,varargin)
            % Constructor of the MovieData object
            %
            % INPUT
            %    channels - a Channel object or an array of Channels
            %    outputDirectory - a string containing the output directory
            %    OPTIONAL - a set of options under the property/key format
            
            if nargin>0
                % Required input fields
                obj.channels_ = channels;
                obj.outputDirectory_ = outputDirectory;
                
                % Construct the Channel object
                nVarargin = numel(varargin);
                if mod(nVarargin,2)==0
                    for i=1 : 2 : nVarargin-1
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
                obj.createTime_ = clock;
            end
        end
        

        %% MovieData specific set/get methods
        function set.movieDataPath_(obj, path)
            % Format the path
            endingFilesepToken = [regexptranslate('escape',filesep) '$'];
            path = regexprep(path,endingFilesepToken,'');
            obj.checkPropertyValue('movieDataPath_',path);
            obj.movieDataPath_=path;
        end
        
        function set.movieDataFileName_(obj, filename)
            obj.checkPropertyValue('movieDataFileName_',filename);
            obj.movieDataFileName_=filename;
        end
        
        function set.channels_(obj, value)
            obj.checkPropertyValue('channels_',value);
            obj.channels_=value;
        end
        
        function set.pixelSize_ (obj, value)
            obj.checkPropertyValue('pixelSize_',value);
            obj.pixelSize_=value;
        end
        
        function set.timeInterval_ (obj, value)
            obj.checkPropertyValue('timeInterval_',value);
            obj.timeInterval_=value;
        end
        
        function set.numAperture_ (obj, value)
            obj.checkPropertyValue('numAperture_',value);
            obj.numAperture_=value;
        end
        
        function set.camBitdepth_ (obj, value)
            obj.checkPropertyValue('camBitdepth_',value);
            obj.camBitdepth_=value;
        end
        
        function set.magnification_ (obj, value)
            obj.checkPropertyValue('magnification_',value);
            obj.magnification_=value;
        end
        function set.binning_ (obj, value)
            obj.checkPropertyValue('binning_',value);
            obj.binning_=value;
        end
        
        function fileNames = getImageFileNames(obj,iChan)
            % Retrieve the names of the images in a specific channel
            
            if nargin < 2 || isempty(iChan), iChan = 1:numel(obj.channels_); end
            assert(all(ismember(iChan,1:numel(obj.channels_))),...
                'Invalid channel numbers! Must be positive integers less than the number of image channels!');
            
            % Delegates the method to the classes
            fileNames = arrayfun(@getImageFileNames,obj.channels_(iChan),...
                'UniformOutput',false);
            if ~all(cellfun(@numel,fileNames) == obj.nFrames_)
                error('Incorrect number of images found in one or more channels!')
            end
        end
        
        function chanPaths = getChannelPaths(obj,iChan)
            %Returns the directories for the selected channels
            if nargin < 2 || isempty(iChan), iChan = 1:numel(obj.channels_); end
            assert(all(ismember(iChan,1:numel(obj.channels_))),...
                'Invalid channel index specified! Cannot return path!');
            
            chanPaths = arrayfun(@(x)obj.channels_(x).channelPath_,iChan,...
                'UniformOutput',false);
        end  
        
        %% Sanitycheck/relocation
        function sanityCheck(obj,varargin)
            % Check the sanity of the MovieData objects
            %
            % First call the superclass sanityCheck. Then call the Channel
            % objects sanityCheck, check image properties and set the 
            % nFrames_ and imSize_ properties. 
            % Save the movie to disk if run successfully
            
            % Call the superclass sanityCheck
            if nargin>1, sanityCheck@MovieObject(obj,varargin{:}); end
            
            % Initialize channels dimensions
            width = zeros(1, length(obj.channels_));
            height = zeros(1, length(obj.channels_));
            nFrames = zeros(1, length(obj.channels_));
            
            % Call subcomponents sanityCheck
            for i = 1: length(obj.channels_)
                [width(i) height(i) nFrames(i)] = obj.channels_(i).sanityCheck(obj);
            end
            
            assert(max(nFrames) == min(nFrames), ...
                'Different number of frames are detected in different channels. Please make sure all channels have same number of frames.')
            assert(max(width)==min(width) && max(height)==min(height), ...
                'Image sizes are inconsistent in different channels.\n\n')
            
            % Define imSize_ and nFrames_;
            if ~isempty(obj.nFrames_)
                assert(obj.nFrames_ == nFrames(1), 'Record shows the number of frames has changed in this movie.')
            else
                obj.nFrames_ = nFrames(1);
            end
            if ~isempty(obj.imSize_)
                assert(obj.imSize_(2) == width(1) && obj.imSize_(1) ==height(1), 'Record shows image size has changed in this movie.')
            else
                obj.imSize_ = [height(1) width(1)];
            end
            
            obj.save();
        end
        
        function relocate(obj,varargin)
            % Relocate the MovieData object
            
            % Run superclass relocate method (for movie path and analysis components)
            [oldRootDir newRootDir]=relocate@MovieObject(obj,varargin{:});
            
            % Relocate the Channel objects
            for i=1:numel(obj.channels_),
                obj.channels_(i).relocate(oldRootDir,newRootDir);
            end
        end
        
        function setFig = edit(obj)
            setFig = movieDataGUI(obj);
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
               status=cellfun(@(x,y) MovieData.checkValue(x,y),property,value);
               return
           end
           
           % Get validator for single property
           validator=MovieData.getPropertyValidator(property);
           propName = regexprep(regexprep(property,'(_\>)',''),'([A-Z])',' ${lower($1)}');
           assert(~isempty(validator),['No validator defined for property ' propName]);
           
           % Return result of validation
           status = isempty(value) || validator(value);
        end
        
        function validator = getPropertyValidator(property) 
            validator = getPropertyValidator@MovieObject(property);
            if ~isempty(validator), return; end
            switch property
                case {'channels_'}
                    validator=@(x) isa(x,'Channel');
                case {'movieDataPath_','movieDataFileName_'}
                    validator=@ischar;
                case {'pixelSize_', 'timeInterval_','numAperture_','magnification_','binning_'}
                    validator=@(x) all(isnumeric(x)) && all(x>0);
                case {'camBitdepth_'}
                    validator=@(x) isposintscalar(x) && ~mod(x, 2);
                otherwise
                    validator=[];
            end
        end
        
        function propName = getPathProperty()
            propName = 'movieDataPath_';
        end
        function propName = getFilenameProperty()
            propName = 'movieDataFileName_';
        end
    end
end