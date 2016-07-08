classdef ExternalNonSingularProcess < ExternalProcess & NonSingularProcess
    %EXTERNALNONSINGULARPROCESS An external process that is also
    %non-singular. This allows for multiple external processes to exist in
    %a MovieObject
    %
    % See also ExternalProcess, NonSingularProcess, Process.run
    
    properties
    end
    
    methods
        function obj = ExternalNonSingularProcess(varargin)
            % See ExternalProcess constructor
            obj = obj@ExternalProcess(varargin{:});
        end
    end
    
end

