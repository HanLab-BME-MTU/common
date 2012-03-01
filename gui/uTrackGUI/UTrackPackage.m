classdef UTrackPackage < TrackingPackage
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
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;

                super_args{1} = owner;              
                super_args{2} = [outputDir filesep 'UTrackPackage'];
                
            end
            % Call the superclass constructor
            obj = obj@TrackingPackage(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'U-Track';
        end

        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = uTrackPackageGUI(varargin{:});
        end
    end
    
end