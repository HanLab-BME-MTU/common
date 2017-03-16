classdef PointSourceDetectionProcess3D < DetectionProcess
%PointSourceDetectionProcess3D is a concrete class of a point source
%detection process for 3d
    
    methods (Access = public)
        function obj = PointSourceDetectionProcess3D(owner, outputDir, funParams)
            % Constructor of the SubResolutionProcess
            super_args{1} = owner;
            super_args{2} = PointSourceDetectionProcess3D.getName;
            super_args{3} = @detectMoviePointSources3D;
            
            if nargin < 3 || isempty(funParams)  % Default funParams
                if nargin <2, outputDir = owner.outputDirectory_; end
                funParams = PointSourceDetectionProcess3D.getDefaultParams(owner,outputDir);
            end
            
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
            obj.is3Dcompatible_ = true;

        end
        
        function movieInfo = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            ip =inputParser;
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output', 'movieInfo',@(x) strcmp(x, 'movieInfo'));
            ip.parse(iChan,varargin{:})
            % z stack!!!!!!
            % Data loading
            s = load(obj.outFilePaths_{1,iChan}, 'movieInfo');
            if ~isempty(s.movieInfo)
                movieInfo = s.movieInfo(ip.Results.iFrame);
            else
                movieInfo = s.movieInfo;
            end
        end
        
        function output = getDrawableOutput(obj)
            output = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Point sources';
            output(1).formatData=@PointSourceDetectionProcess3D.formatOutput;
        end
        
        function drawImaris(obj,iceConn)
            
            dataSet = iceConn.mImarisApplication.GetDataSet;
            
            nChan = numel(obj.owner_.channels_);
            
            %Create data container
            spotFolder = iceConn.mImarisApplication.GetFactory.CreateDataContainer;
            spotFolder.SetName(obj.name_);
            iceConn.mImarisApplication.GetSurpassScene.AddChild(spotFolder,-1);
            for iChan = 1:nChan
                
                
                if obj.checkChannelOutput(iChan)
                    %TEMP - doesn't support timelapse yet!!
                    pts = obj.loadChannelOutput(iChan);
                    if ~isempty(pts)
                        chanCol = iceConn.mapRgbaScalarToVector(dataSet.GetChannelColorRGBA(iChan-1));
                        nPts = numel(pts.x);
                        ptXYZ = [pts.y' pts.x' pts.z'];%Swap xy for imaris display
                        ptXYZ(:,1:2) =ptXYZ(:,1:2) * obj.owner_.pixelSize_ / 1e3;
                        ptXYZ(:,3) =ptXYZ(:,3) * obj.owner_.pixelSizeZ_ / 1e3;
                        ptRad = pts.s(1,:)' * obj.owner_.pixelSize_ / 1e3;
                        ptObj = iceConn.createAndSetSpots(ptXYZ,zeros(nPts,1),...
                            ptRad,['Spots ' char(dataSet.GetChannelName(iChan-1))],chanCol,spotFolder);                    

                        if ismethod(ptObj,'SetRadiiXYZ')
                            %Make sure we have a version of imaris which supports anisotropic points                        
                            ptRad = zeros(nPts,3);
                            ptRad(:,1:2) = repmat(pts.s(1,:)' * obj.owner_.pixelSize_ / 1e3,[1 2]);
                            ptRad(:,3) = pts.s(2,:) * obj.owner_.pixelSizeZ_ / 1e3;
                            ptObj.SetRadiiXYZ(ptRad);                                            
                        end
                    end
                end
            end            
            
        end

        %% TODO WIP -- TO FINISH (using Philippe'method)
        function scales = getEstSigmaPSF3D(obj, channel, varargin)

            volList=[];
            for i=1:5
                volList = [volList double(obj.owner_.getChannel(channel).loadStack(i))];
            end
            scales = getGaussianPSFsigmaFrom3DData(volList, varargin{:});
            disp(['Estimated scales: ' num2str(scales)]);

        end

    end
    methods (Static)
        
        function name = getName()
            name = 'Point source detection 3D';
        end
        function h = GUI()
            h = @pointSourceDetectionProcessGUI3D;
        end
        
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1;
            
            %## TODO
            funParams.InputImageProcessIndex = 0; % ?? (can we add some way to check what is availble.)
            
            funParams.MaskChannelIndex = []; %1:numel(owner.channels_);
            funParams.MaskProcessIndex = [];            
            funParams.OutputDirectory = [outputDir  filesep 'point_sources'];
            
            funParams.alpha=.05;
            funParams.Mode = {'xyzAc'};
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RemoveRedundant = true;
            funParams.RedundancyRadius = .25;
            funParams.RefineMaskLoG = false;
            funParams.RefineMaskValid = false;
            % funParams.ConfRadius = []; 
            % funParams.WindowSize = [];       

            % sigma estimation            
            nChan = numel(owner.channels_);
            funParams.filterSigma = 1.2*ones(1,nChan); %Minimum numerically stable sigma is ~1.2 pixels.
            hasPSFSigma = arrayfun(@(x) ~isempty(x.psfSigma_), owner.channels_);
            funParams.filterSigma(hasPSFSigma) = [owner.channels_(hasPSFSigma).psfSigma_];            
            funParams.filterSigma(funParams.filterSigma<1.2) = 1.2; %Make sure default isn't set to too small.
            funParams.filterSigma = repmat(funParams.filterSigma,[2 1]); %TEMP - use z-PSF estimate as well!!
            % For the GUI
            funParams.filterSigmaXY = funParams.filterSigma(1,1);
            funParams.filterSigmaZ = funParams.filterSigma(2,1);
            funParams.ConfRadius = arrayfun(@(x)(2*x), funParams.filterSigma);
            funParams.WindowSize = arrayfun(@(x)(ceil(4*x)), funParams.filterSigma);                       

            %list of parameters which can be specified at a per-channel
            %level. If specified as scalar these will  be replicated
            funParams.PerChannelParams = {'alpha','Mode','FitMixtures','MaxMixtures','RedundancyRadius',...
                'ConfRadius','WindowSize','RefineMaskLoG','filterSigma'};
            funParams = prepPerChannelParams(funParams, nChan);
        end
        
        function positions = formatOutput(pstruct)
            positions = formatOutput@DetectionProcess(pstruct);
            %positions = positions(pstruct.isPSF, :);
        end
          

    end    
end