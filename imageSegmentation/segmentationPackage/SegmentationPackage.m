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
                super_args{2} = SegmentationPackage.getName;
                % Dependency Matrix (same length as process class name
                % string)
                super_args{3} = SegmentationPackage.getDependencyMatrix;
                
                % Process CLASS NAME string (same length as dependency matrix)
                % Must be accurate process class name
                segmentationClasses = {
                    @SegmentationProcess,...
                    @MaskRefinementProcess};
                super_args{4} = cellfun(@func2str,segmentationClasses,...
                    'UniformOutput',false);
                
                super_args{5} = [outputDir filesep 'SegmentationPackage'];
                
            end
            % Call the superclass constructor
            obj = obj@Package(super_args{:},'processClassHandles_',segmentationClasses);
        end

        function processExceptions = sanityCheck(obj,varargin)
            
            % Check that the channels have a psf function
            nProc = length(obj.processClassNames_);
            ip = inputParser;
            ip.CaseSensitive = false;
            ip.addRequired('obj');
            ip.addOptional('full',true, @(x) islogical(x));
            ip.addOptional('procID',1:nProc,@(x) all(ismember(x,1:nProc)) || strcmp(x,'all'));
            ip.parse(obj,varargin{:});
            full = ip.Results.full;
            procID = ip.Results.procID;
            if strcmp(procID,'all'), procID = 1:nProc;end
            
            if full
                validProc = procID(~cellfun(@isempty,obj.processes_(procID)));
                if all(ismember([1 2],validProc))
                    % Find the segmentation process index and set it in the
                    % mask refinement process
                    funParams.SegProcessIndex = find(cellfun(@(x) isequal(x,obj.processes_{1}),...
                        obj.owner_.processes_));
                    parseProcessParams(obj.processes_{2},funParams);
                end
            end
            processExceptions = sanityCheck@Package(obj,varargin{:});
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
        
        function name = getName()
            name = 'Segmentation';
        end
        
        function varargout = start(varargin)
            % Start the package GUI
            varargout{1} = segmentationPackageGUI(varargin{:});
        end
    end
    
end

