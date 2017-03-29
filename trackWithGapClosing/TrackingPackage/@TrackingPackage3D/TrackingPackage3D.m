classdef TrackingPackage3D < Package
    % An abstract class for a geeneric Tracking Package
    
    methods
        function obj = TrackingPackage3D(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'TrackingPackage3D'];
            end
                 
            % Call the superclass constructor
            obj = obj@Package(super_args{:});        
        end
        
        
        function [status, processExceptions] = sanityCheck(obj, varargin) % throws Exception Cell Array

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


            % Check channels?
            % Check detections process when feeding to tracking ...

        end    
        
    end

    methods (Static)
        function name = getName()
            name = 'U-Track 3D';
        end 

        function m = getDependencyMatrix(i,j)   
                %1 2 3 4 5 6 7 8
            m = [0 0 0 0 0 0 0 0;   %1 MIP Process (1) [optional]
                 0 0 0 0 0 0 0 0;   %2 DetectionProcess
                 0 1 0 0 0 0 0 0;   %3 Registration [optional] (if not done lab ref is selected)
                 0 0 0 0 0 0 0 0;   %4 Fiduciary/Pole Detection Process [optional?]
                 0 1 2 2 0 0 0 0;   %5 TrackingProcess (can use spindle, stage, or lab ref)
                 0 1 2 2 1 0 0 0;   %6 PostTrackingProcess
                 0 1 2 2 1 0 0 0;   %7 Projection Process [optional]
                 0 1 2 2 1 2 2 0;]; %8 export external rendering Process [optional]?

            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI 
            varargout{1} = TrackingPackage3DGUI(varargin{:});
        end

        function classes = getProcessClassNames(index)
            classes = {'ExternalProcess'    ,...
                       'DetectionProcess'   ,...
                       'StageDriftCorrectionProcess'    ,...
                       'ExternalProcess'    ,...
                       'TrackingProcess'    ,...
                       'PostTrackingProcess',...
                       'ExternalProcess'    ,...
                       'ExternalProcess'     ...
                       };
            if nargin==0, index=1:numel(classes); end
            classes=classes(index);
        end
        
        function objects = getConcretePackages(varargin)
            % If input, check if 2D or 3D movie(s).
            ip = inputParser;
            ip.addOptional('MO', [], @(x) isa(x,'MovieData') || isa(x,'MovieList'));
            ip.parse(varargin{:});
            MO = ip.Results.MO;
            
            if ~isempty(MO)
                if isa(MO,'MovieList')
                    MD = MO.getMovie(1);
                elseif length(MO) > 1
                    MD = MO(1);
                else
                    MD = MO;
                end                
            end

            % Note this will redirect to superclass TrackingPackage, Not this superclass [TrackingPackage3D]
            % if 2D movie detected...
            % 2D options
            objects(1).name = '[2D]Single particles';
            objects(1).packageConstr = @UTrackPackage;
            objects(2).name = '[2D]Microtubules plus-ends';
            objects(2).packageConstr = @PlusTipTrackerPackage;                        
            objects(3).name = '[2D]Nuclei';
            objects(3).packageConstr = @NucleiTrackingPackage;
            % 3D options
            %% NOTE: for coding expediency making new TrackingPackage3D 
            %% Need to make sure this will not cause errors...
            objects(4).name = '[3D]Single particles';
            objects(4).packageConstr = @UTrackPackage3D; 
            
            % Note this will redirect to superclass TrackingPackage3D, Not this superclass [TrackingPackage]
            %% Need to make sure this will not cause errors...
            objects(5).name = '[3D]Microtubules plus-ends';
            objects(5).packageConstr = @PlusTipTrackerPackage3D;             

            
            if isempty(MD)
               % disp('MovieData properties not specified (2D vs. 3D)...showing only 3D options');
               % Assume user knowlingly called 3D package.
                objects(1:3) = [];
            elseif MD.is3D
                objects(1:3) = [];
            elseif ~MD.is3D
                warning('Detected 2D movie');
                warning('Redirecting to 2D Package Constructors (TrackingPackage superclass)');
                objects(4:5) = [];
            end
        end
        
    end 
end