classdef UTrackPackage3D < TrackingPackage
    % A concrete process for UTrack Package
    
    methods (Access = public)
        function obj = UTrackPackage3D(varargin)

            % Call the superclass constructor
            obj = obj@TrackingPackage(varargin{:});
        end        

        function [status processExceptions] = sanityCheck(obj, varargin) % throws Exception Cell Array
            
            %% TODO - add more to sanitycheck
            disp('TODO: SanityCheck!');
            missingMetadataMsg = ['Missing %s! The %s is necessary to analyze '...
            '3D Tracking Movies. Please edit the movie and fill the %s.'];
            errorMsg = @(x) sprintf(missingMetadataMsg, x, x, x);
            
            assert(obj.owner_.is3D, errorMsg('MovieData is not 3D!'));
            assert(~isempty(obj.owner_.pixelSize_), errorMsg('pixel size not defined!'));
            assert(~isempty(obj.owner_.pixelSizeZ_), errorMsg('pixel Z size defined!'));
            assert(~isempty(obj.owner_.timeInterval_), errorMsg('time interval defined!'));
            [status, processExceptions] = sanityCheck@Package(obj, varargin{:});

            % possible PSF sanity check?

        end
    end

    methods (Static)

        function name = getName()
            name = 'UTrackPackage3D';
        end

        function varargout = GUI(varargin)
            varargout{1} = uTrackPackageGUI(varargin{:});
        end

        function procConstr = getDefaultProcessConstructors(index)
            procConstr = {
                @(x,y)PointSourceDetectionProcess3D(x,y,UTrackPackage3D.getDefaultDetectionParams(x,y)),...
                @(x,y)TrackingProcess(x,y,UTrackPackage.getDefaultTrackingParams(x,y))};%,...

                % @MotionAnalysisProcess};
            if nargin==0, index=1:numel(procConstr); end
            procConstr=procConstr(index);
        end
        
        function funParams = getDefaultDetectionParams(owner, outputDir)

            funParams = PointSourceDetectionProcess3D.getDefaultParams(owner, outputDir);
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'pointsource3D_detect'];

            %% TODO -- Update 
            iProc = owner.getProcessIndex('MaskRefinementProcess', 'askUser', false);
            if isempty(iProc)
                disp('Note: No Cell Segmentation Mask found');
                funParams.MaskChannelIndex = []; %1:numel(owner.channels_);
                funParams.MaskProcessIndex = [];            
            else
                funParams.MaskProcessIndex = iProc; % Specify Process Index with cell mask output
                funParams.MaskChannelIndex = 1:numel(owner.channels_);
            end
            
            funParams.alpha = .01;

            % Extra outputs? (convert to obj methods)
            % funParams.printAll ???? <<<<<<<<<< (convert to obj methods)
            % funParams.showAll ???? <<<<<<< (convert to obj methods)

            funParams = prepPerChannelParams(funParams, numel(owner.channels_));
        end

        function funParams = getDefaultTrackingParams(owner,outputDir)
            funParams = TrackingProcess.getDefaultParams(owner,outputDir);

            % Set default kalman functions
            funParams.kalmanFunctions = TrackingProcess.getKalmanFunctions(1);

            % Set default cost matrices
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
        end
        
    end
    
end