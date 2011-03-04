classdef SubResolutionProcess < DetectionProcess
    %
    % Chuangang Ren
    % 11/2010

    properties(SetAccess = protected, GetAccess = public)

        filenameBase_ % cell array storing the base of file name 
        digits4Enum_ % cell array stroring the number of digits for frame enumeration (1-4)
        filename_ % file name of result data
    end    
    
    methods (Access = public)
        function obj = SubResolutionProcess(owner, outputDir, channelIndex, funParams )

            super_args{1} = owner;
            super_args{2} = 'Sub-Resolution';
            super_args{3} = @detectSubResFeatures2D_StandAlone;
            
            if nargin < 2 || isempty(outputDir)
                outputDir = owner.outputDirectory_ ; % package folder
            end                
            
            if nargin < 3 || isempty(channelIndex) % Default channel Index
                channelIndex = 1:length(owner.channels_);
            end
            
            super_args{4} = channelIndex;

            
            if nargin < 4 || isempty(funParams)  % Default funParams                         
                    
                % movieParam
                funParams.movieParam.imageDir = owner.channels_(channelIndex(1)).channelPath_; % Note: channel-specific
                funParams.movieParam.filenameBase = []; % Note: channel-specific
                funParams.movieParam.firstImageNum = 1;
                funParams.movieParam.lastImageNum = owner.nFrames_;
                funParams.movieParam.digits4Enum = []; % Note: channel-specific
                
                % detectionParam
%                 funParams.detectionParam.psfSigma = [];
%                 funParams.detectionParam.bitDepth = owner.camBitdepth_;
                funParams.detectionParam.alphaLocMax = .05;
                funParams.detectionParam.integWindow = 0;
                funParams.detectionParam.doMMF = 0;
                funParams.detectionParam.testAlpha = struct('alphaR', .05,'alphaA', .05, 'alphaD', .05,'alphaF',0);
                funParams.detectionParam.numSigmaIter = 0;
                funParams.detectionParam.visual = 0;
                funParams.detectionParam.background = []; 
                
                % saveResults
                funParams.saveResults.dir = [outputDir  filesep 'Sub_Resolution_Detection' filesep];
                funParams.saveResults.filename = []; % Note: channel-specific
                
                % Set up psfSigma and bitDepth
                na = owner.numAperature_;
                ps = owner.pixelSize_;
                wl = owner.channels_(1).emissionWavelength_;
                bd = owner.camBitdepth_;
                
                if ~isempty( na ) && ~isempty( ps ) && ~isempty( wl )
                    funParams.detectionParam.psfSigma = 0.21*wl/na/ps;
                else
                    funParams.detectionParam.psfSigma = [];
                end
                
                if ~isempty(bd)
                    funParams.detectionParam.bitDepth = bd;
                else
                    funParams.detectionParam.bitDepth = [];
                end
                
            end

            super_args{5} = funParams;

            obj = obj@DetectionProcess(super_args{:});    
            
            % Visual parameters ( Default: channel 1 )
            obj.visualParams_.startend = [1 owner.nFrames_];
            obj.visualParams_.saveMovie = 1;
            obj.visualParams_.movieName = [];
            obj.visualParams_.dir2saveMovie = funParams.saveResults.dir;            
            obj.visualParams_.filterSigma = 0;
            obj.visualParams_.showRaw = 1;
            obj.visualParams_.intensityScale = 1;
            file = owner.getImageFileNames(1);
            obj.visualParams_.firstImageFile = [owner.channels_(1).channelPath_ filesep file{1}{1}];

            
            % Get file name base and digits for enumeration
            [obj.filenameBase_ obj.digits4Enum_] = SubResolutionProcess.getFilenameBody(owner);
            obj.filename_ = 'detection_result.mat';
        end    
        
        % Set result file name
        function setFileName(obj, name)
           obj.filename_ = name; 
        end
        
        function runProcess(obj)
        % Run the process!
            for i = obj.channelIndex_
                
                obj.funParams_.movieParam.imageDir = [obj.owner_.channels_(i).channelPath_ filesep];
                obj.funParams_.movieParam.filenameBase = obj.filenameBase_{i};
                obj.funParams_.movieParam.digits4Enum = obj.digits4Enum_{i};
                obj.funParams_.saveResults.filename = ['Channel_' num2str(i) '_' obj.filename_];
                
                %Check/create directory
                if ~exist(obj.funParams_.saveResults.dir,'dir')
                    mkdir(obj.funParams_.saveResults.dir)
                end
                
                if ~obj.overwrite_
                    % file name enumaration
                    obj.funParams_.saveResults.filename = enumFileName(obj.funParams_.saveResults.dir, obj.funParams_.saveResults.filename);
                end
                
                % Test (commentable)
%                 obj.funParams_.movieParam,obj.funParams_.detectionParam,obj.funParams_.saveResults

                [obj.outParams_{i}.movieInfo, obj.outParams_{i}.exceptions, obj.outParams_{i}.localMaxima, ...
                    obj.outParams_{i}.background, obj.outParams_{i}.psfSigma] = ...
                    obj.funName_(obj.funParams_.movieParam, obj.funParams_.detectionParam, obj.funParams_.saveResults );
            end
        end
        
    end
    methods (Static)
        
        function [base digits4Enum]= getFilenameBody(owner)
            % Get the base of file name for all channels in movie data
            % "owner"
            
            fileNames = owner.getImageFileNames;
            base = cell(1, length(owner.channels_));
            digits4Enum = cell(1, length(owner.channels_));
            
            for i = 1 : length(owner.channels_)
                
                [x1 base{i} digits4Enum{i} x4] = getFilenameBody(fileNames{i}{1});
                digits4Enum{i} = length(digits4Enum{i});
            end
            
        end    
        
    end
    
end