classdef UTrackPackage < Package
    % A concrete process for UTrack Package
    
    methods (Access = public)
        function obj = UTrackPackage (owner,varargin)
            % Construntor of class MaskProcess
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieObject'));
                ip.addOptional('outputDir',owner.outputDirectory_,@(x) isa(x,'MovieObject'));
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;

                super_args{1} = owner;              
                super_args{2} = [outputDir filesep 'UTrackPackage'];
                
            end
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'U-Track';
        end
        
        function m = getDependencyMatrix(i,j)

            m = [0 0; %1 DetectionProcess
                1 0]; %2 TrackingProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end

        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = uTrackPackageGUI(varargin{:});
        end
               
        function procConstr = getDefaultProcessConstructors(index)
            uTrackProcConstr = {
                @SubResolutionProcess,...
                @TrackingProcess};
            
            if nargin==0, index=1:numel(uTrackProcConstr); end
            procConstr=uTrackProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            uTrackClasses = {
                'DetectionProcess',...
                'TrackingProcess'};
            if nargin==0, index=1:numel(uTrackClasses); end
            classes=uTrackClasses(index);
        end
        
    end
    
end