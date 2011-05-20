classdef SegmentationPackage < Package
    % A concrete process for Segmentation Package
    
    methods (Access = public)
        function obj = SegmentationPackage (owner,outputDir)
            % Construntor of class MaskProcess
            if nargin == 0
                super_args = {};
            else
                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = 'Segmentation';
                % Dependency Matrix (same length as process class name
                % string)
                super_args{3} = SegmentationPackage.getDependencyMatrix;
                
                % Process CLASS NAME string (same length as dependency matrix)
                % Must be accurate process class name
                segmentationClasses = {
                    @ThresholdProcess,...
                    @MaskRefinementProcess};
                super_args{4} = cellfun(@func2str,segmentationClasses,...
                    'UniformOutput',false);
                
                super_args{5} = [outputDir filesep 'SegmentationPackage'];
                
            end
            % Call the superclass constructor
            obj = obj@Package(super_args{:},'processClassHandles_',segmentationClasses);
        end
        
    end
    methods (Static)
        
        function m = getDependencyMatrix()
            % Get dependency matrix
            m = [0 0; % SegmentationProcess
                1 0]; % MaskRefinementProcess
        end
        
        function id = getOptionalProcessId()
            % Get the optional process id
            id = [];
        end
        
        function varargout = start(varargin)
            % Start the package GUI
            varargout{1} = segmentationPackageGUI(varargin{:});
        end
    end
    
end

