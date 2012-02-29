classdef PostTrackingProcess < DataProcessingProcess
    % An abstract class associated to the tracks post-tracking (e.g. classification)
    %
    % Sebastien Besson, Feb 2012
    
    methods (Access = public)
        function obj = PostTrackingProcess(owner, varargin)
            obj = obj@DataProcessingProcess(varargin{:});
        end
    end
    methods (Static)
        function procClasses = getConcreteClasses()
            procClasses = ...
                {'CometPostTrackingProcess';};
        end

    end    
end