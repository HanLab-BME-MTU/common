classdef UTrackPackage3D < TrackingPackage
    % A concrete process for UTrack Package
    
    methods (Access = public)
        function obj = UTrackPackage3D (varargin)
            % Call the superclass constructor
            obj = obj@TrackingPackage(varargin{:});
        end        

        function [status processExceptions] = sanityCheck(obj, varargin) % throws Exception Cell Array
            %% TODO - sanitycheck
            disp('TODO: SanityCheck!');

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
                @(x,y)PointSourceDetectionProcess(x,y,UTrackPackage3D.getDefaultDetectionParams(x,y)),...
                @(x,y)TrackingProcess(x,y,UTrackPackage.getDefaultTrackingParams(x,y)),...
                @MotionAnalysisProcess};
            if nargin==0, index=1:numel(procConstr); end
            procConstr=procConstr(index);
        end
        
        function funParams = getDefaultDetectionParams(owner, outputDir)

            funParams = PointSourceDetectionProcess.getDefaultParams(owner, outputDir);
            
            %% TODO - Verify ideal default settings here.
            % Set default parameters
            funParams.ChannelIndex = 1;

            % Check if segmentation occured.
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

            funParams.OutputDirectory = [outputDir  filesep 'point_sources'];
            funParams.alpha=.05;
            funParams.maskRadius=40;
            funParams.Mode = {'xyAc'};
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RedundancyRadius = .25;
            funParams.UseIntersection = true;            
            funParams.PreFilter = true;
            %list of parameters which can be specified at a per-channel
            %level. If specified as scalar these will  be replicated
            funParams.PerChannelParams = {'alpha','Mode','FitMixtures','MaxMixtures','RedundancyRadius','filterSigma','PreFilter','ConfRadius','WindowSize'};
            
            nChan = numel(owner.channels_);
            funParams.filterSigma = 1.2*ones(1,nChan); %Minimum numerically stable sigma is ~1.2 pixels.
            hasPSFSigma = arrayfun(@(x) ~isempty(x.psfSigma_), owner.channels_);
            funParams.filterSigma(hasPSFSigma) = [owner.channels_(hasPSFSigma).psfSigma_];            
            funParams.filterSigma(funParams.filterSigma<1.2) = 1.2; %Make sure default isn't set to too small.
            
            funParams.ConfRadius = arrayfun(@(x)(2*x),funParams.filterSigma);
            funParams.WindowSize = arrayfun(@(x)(ceil(4*x)),funParams.filterSigma);
            
            funParams = prepPerChannelParams(funParams,nChan);
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