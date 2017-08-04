classdef  ComputeMIPProcess < ImageProcessingProcess
    % Concrete class for a computing Maximum Intensity Projections (MIP)
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

                super_args{1} = owner;
                super_args{2} = ComputeMIPProcess.getName;
                super_args{3} = @computeMovieMIP;
                if isempty(funParams)
                    funParams = ComputeMIPProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@ImageProcessingProcess(super_args{:});
            obj.is3Dcompatible_ = false; % outputs are 2D
        end

        function h = draw(obj, varargin)
            % Function to draw process output
            outputList = obj.getDrawableOutput();  
                               
            ip = inputParser;
            ip.addRequired('obj',@(x) isa(x,'Process'));
            ip.addRequired('iChan',@isnumeric);
            ip.addOptional('iFrame',[],@isnumeric);
            ip.addOptional('iZ',[], @(x) ismember(x,1:obj.owner_.zSize_));
            ip.addParameter('output', [], @(x) all(ismember(x,{outputList.var})));
            ip.KeepUnmatched = true;
            ip.parse(obj, varargin{:});

            iOutput = find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if iOutput == 2
                iOutput_XY = 1;
            end
            % Display all channels
            switch ip.Results.output
                case 'merged'
                    if numel(obj.owner_.channels_) > 1, cdim=3; else cdim=1; end
                        data = zeros([obj.owner_.imSize_ cdim]);

                    for iChan = 1:numel(obj.owner_.channels_)
                        imData = obj.loadChannelOutput(iChan, ip.Results.iFrame, iOutput);
                        data(:,:,iChan) = outputList(iOutput).formatData(imData);
                    end                  

                    try
                        assert(~isempty(obj.displayMethod_{iOutput,1}));
                    catch ME
                        obj.displayMethod_{iOutput,1}=...
                            outputList(iOutput).defaultDisplayMethod();
                    end

                    % Create graphic tag and delegate drawing to the display class
                    tag = ['process' num2str(obj.getIndex()) '_MergedOutput'];
                    h = obj.displayMethod_{iOutput}.draw(data, tag, ip.Unmatched);

                case {'XY'}
                    % Call superclass method
                    % imData = obj.loadChannelOutput(iChan, ip.Results.iFrame, 'iOutput', iOutput_XY);
                    % data(:,:,iChan) = outputList(iOutput).formatData(imData);
                    h = draw@ImageProcessingProcess(obj, [varargin, 'iOutput', iOutput_XY]);
                
                case {'ZY', 'ZX'}

                    % h = draw@ImageProcessingProcess(obj, [varargin, 'iOutput', iOutput]);
                    h = draw@ImageProcessingProcess(obj, varargin{1}, varargin{2},...
                            varargin{3:end}, 'iOutput', iOutput);
                    %     h = draw@ImageProcessingProcess(obj, varargin{:}, 'iOutput', iOutput);
                
                % case 'ZX'
                
                %     h = draw@ImageProcessingProcess(obj, varargin{:}, 'iOutput', iOutput);
                
                otherwise
                    error('Incorrect Output Var type');
            end
        end
        
        function output = getDrawableOutput(obj, varargin)
            
            output(1).name = 'Merged';
            output(1).var = 'merged';
            output(1).formatData = @mat2gray;
            output(1).defaultDisplayMethod = @ImageDisplay;
            output(1).type = 'image';
            
            n = length(output)+1;
            output(n).name = 'XY';
            output(n).var = 'XY';
            output(n).formatData = @mat2gray;
            output(n).defaultDisplayMethod = @ImageDisplay;
            output(n).type = 'image';
            
            n = length(output)+1;
            output(n).name = 'ZY';
            output(n).var = 'ZY';
            output(n).formatData = @mat2gray;
            output(n).defaultDisplayMethod = @ImageDisplay;
            output(n).type = 'image';
            
            n = length(output)+1;
            output(n).name = 'ZX';
            output(n).var = 'ZX';
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
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'MIP'];
        end
    end
end