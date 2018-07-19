classdef IntegratorPackage < Package
    % The main class of the Integrator package
    
    % Sebastien Besson, Sep 2011
    
    methods
        function obj = IntegratorPackage(owner,varargin)
            % Constructor of class QFSMPackage
            
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieList'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;

                % Owner: MovieData object
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'IntegratorPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        function [status processExceptions] = sanityCheck(obj,varargin)

            movies = obj.owner_.getMovies;
            nMovies=numel(movies);
            
            % Check that the frame rate is input for all movies
            assert(all(cellfun(@(x) ~isempty(x.timeInterval_),movies)),...
                'Please fill the time interval for all movies.');
            
            % Check signal input consistency in-between movies
            signalInput=cell(nMovies,1);
            for i=1:nMovies,
                try 
                    signalInput{i}=movies{i}.getSampledOutput;
                catch ME
                    throw(MException('lccb:ml:sanitycheck','Movie %g\n\n%s',...
                        i, ME.message));
                end
            end
            nInput=unique(cellfun(@numel,signalInput));
            assert(isscalar(nInput),'Different number of sampled signal per movie');
            for i=1:nInput, 
                signalName = unique(cellfun(@(x) x(i).name,signalInput,'Unif',false));
                assert(numel(signalName)==1,'Different signal type in betwen movies');
            end
            
            % Call superclass sanityCheck
            [status processExceptions] = sanityCheck@Package(obj,varargin{:});
        end
        
    end
    
    methods (Static)
        
        function m = getDependencyMatrix(i,j)
            
            m = [0 0 0;  %1 SignalPreprocessingProcess
                1 0 0;]; %2 SignalProcessingProcess

            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function name = getName()
            name='Image BLAST';
        end

        function name = getMovieClass()
            name='MovieList';
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = integratorPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            integratorProcConstr = {
                @SignalPreprocessingProcess,...
                @SignalProcessingProcess};
            
            if nargin==0, index=1:numel(integratorProcConstr); end
            procConstr=integratorProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            integratorClasses = {
                'SignalPreprocessingProcess',...
                'SignalProcessingProcess'};
            if nargin==0, index=1:numel(integratorClasses); end
            classes=integratorClasses(index);
        end
    end
end
