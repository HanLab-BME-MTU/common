classdef UTrackPackage < Package
% A concrete process for UTrack Package
    
methods (Access = public)
    function obj = UTrackPackage (owner,outputDir)
           % Construntor of class MaskProcess
           if nargin == 0
              super_args = {};
           else
               % Owner: MovieData object
               super_args{1} = owner;
               super_args{2} = 'U-Track'; 
               % Dependency Matrix (same length as process class name
               % string)
               super_args{3} = UTrackPackage.getDependencyMatrix;
                                
               % Process CLASS NAME string (same length as dependency matrix)
               % Must be accurate process class name
               uTrackClasses = {
                   @DetectionProcess,...
                   @TrackingProcess};
               super_args{4} = cellfun(@func2str,uTrackClasses,...
                   'UniformOutput',false);
                            
               super_args{5} = [outputDir filesep 'UTrackPackage'];
                
           end
           % Call the superclass constructor 
           obj = obj@Package(super_args{:},'processClassHandles_',uTrackClasses);
    end 

end
methods (Static)
    
        function m = getDependencyMatrix()
            % Get dependency matrix
               m = [0 0;
                    1 0];
        end        
        
        function id = getOptionalProcessId()
            % Get the optional process id
            id = [];
        end
        
        function varargout = start(varargin)
            % Start the package GUI
            varargout{1} = uTrackPackageGUI(varargin{:});
        end
        
end
    
end