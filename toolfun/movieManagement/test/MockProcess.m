classdef MockProcess < Process
    methods (Access = public)
        function obj = MockProcess(owner, varargin)
            obj = obj@Process(owner, MockProcess.getName);
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'Mock process';
        end
        
        function funParams = getDefaultParams(owner)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.parse(owner)
            
            % Set default parameters
            funParams.MockParam1 = true;
        end
    end
end