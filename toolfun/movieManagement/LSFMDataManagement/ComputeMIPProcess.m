classdef  ComputeMIPProcess < ImageProcessingProcess
    % Concrete class for a computing Maximum Intensity Projections (MIP)
    %
    % Andrew R. Jamieson Aug. 2017
    
    methods
        function obj = ComputeMIPProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Constructor for EfficientSubpixelRegistrationProcess

                super_args{1} = owner;
                super_args{2} = ComputeMIPProcess.getName;
                super_args{3} = @computeMovieMIP;
                if isempty(funParams)
                    funParams = ComputeMIPProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@ImageProcessingProcess(super_args{:});
        end
        
        function output = getDrawableOutput(obj, varargin)
            output = getDrawableOutput@ImageProcessingProcess(obj);
            n = length(output)+1;
            output(n).name = 'Merged';
            output(n).var = 'merged';
            output(n).formatData = @mat2gray;
            output(n).defaultDisplayMethod = @ImageDisplay;
            output(n).type = 'image';
        end
    end
    
    methods (Static)
        function name = getName()
            name = 'Maximum Intensity Projection';
        end

        % function h = GUI()
        %     % h = @ComputeMIPProcessGUI;
        %     h = @abstractProcessGUI;
        % end
        
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'MIP'];
        end
    end
end