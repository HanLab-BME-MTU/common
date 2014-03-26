classdef ExampleProcess < Process
    methods (Access = public)
        function obj = ExampleProcess(owner)
            obj = obj@Process(owner, ExampleProcess.getName);
            obj.funName_ = @wrapperMovieExample;
            obj.funParams_ = ExampleProcess.getDefaultParams(owner);
        end
        
        
        function output = loadChannelOutput(obj, iChan)
            s = load(obj.outFilePaths_{iChan});
            output = s.Imean;          
        end
    end
    methods (Static)
        function name = getName()
            name = 'Example';
        end
        
        function funParams = getDefaultParams(owner)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieObject'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner)
            
            % Set default parameters
            funParams.OutputDirectory = [ip.Results.outputDir  filesep 'stats'];        end
    end
end