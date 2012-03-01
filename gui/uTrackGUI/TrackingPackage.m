classdef TrackingPackage < Package
    % A concrete process for Biosensor Package
    methods
        function obj = TrackingPackage(varargin)
            % Call the superclass constructor
            obj = obj@Package(varargin{:});        
        end
    end
    methods (Static)
        function m = getDependencyMatrix(i,j)   
            m = [0 0 0;  %1 DetectionProcess
                 1 0 0;  %2 TrackingProcess
                 0 1 0;];%3 MotionAnalysisProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = trackingPackageGUI(varargin{:});
        end
        function classes = getProcessClassNames(index)
            classes = {
                'DetectionProcess',...
                'TrackingProcess',...
                'MotionAnalysisProcess'};
            if nargin==0, index=1:numel(classes); end
            classes=classes(index);
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            procConstr = {
                @SubResolutionProcess,...
                @TrackingProcess,...
                @MotionAnalysisProcess};
            if nargin==0, index=1:numel(procConstr); end
            procConstr=procConstr(index);
        end
        
    end 
end