classdef SegmentationPackage < Package
% A concrete process for BioSensor Package
    
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
               super_args{4} = {'SegmentationProcess'};
                            
               super_args{5} = [outputDir filesep 'SegmentationPackage'];
                
           end
           % Call the superclass constructor 
           obj = obj@Package(super_args{:});
    end
    
    function processExceptions = sanityCheck(obj,full,procID) % throws Exception Cell Array
        % Sanity Check
        % full package panity check: true or false
        nProcesses = length(obj.processClassNames_);
            
        if nargin < 2
            full = true;
            procID = 1:nProcesses;
        end
            
        if nargin < 3
           procID = 1:nProcesses ;
        end
            
        if strcmp(procID,'all')
            procID = 1:nProcesses;
        end
            
        if any(procID > nProcesses)
            error('User-defined: process id exceeds number of processes');
        end
            
        processExceptions = obj.checkProcesses(full,procID);  % throws Exception Cell Array
                     
    end
    

end
methods (Static)
        
        function m = getDependencyMatrix()
            % Get dependency matrix
               m = 0;
        end        
        
        function id = getOptionalProcessId()
            % Get the optional process id
            id = [];
        end
end
    
end

