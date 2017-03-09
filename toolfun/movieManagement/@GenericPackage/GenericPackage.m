classdef GenericPackage < Package
    %GenericPackage Creates a Generic package out of already constructed
    %processes
    %
    % INPUT
    % GenericPackage(MovieObject) - Creates a generic package out of the
    %       MovieObject's current processes
    % GenericPackage(Processes) - Creates a generic package out of certain
    %       processes in a cell array. Assumes owner is
    %       processes{1}.getOwner()
    % GenericPackage( __ , outputDirectory) - Sets the outputDirectory_
    %       property, effect unclear
    % GenericPackage( __ , 'name_', string) - Sets the name of the
    %       GenericPackage
    % GenericPackage( __ , 'dependencyMatrix_', matrix) - Sets the
    %       dependency matrix. Default: diag(ones(length(obj.processes_)-1,1),-1)
    %
    %
    % USAGE
    % pkg = GenericPackage(MD,MD.outputDirectory_,'name_','Hello World Package')
    % pkg.dependencyMatrix_ = tril(ones(length(pkg.processes_)),-1);
    % MD.addPackage(pkg);
    % % Invoke GUI
    % MD.packages_{1}.GUI(MD);
    %
    % See also ExternalProcess, cliGUI
    
    % Note: This exploits the fact that Static Abstract methods do not have
    % be implemented as Static
     
    % Mark Kittisopikul, March 2017
    % Jaqaman Lab
    % UT Southwestern
    
    properties
        name_ = 'GenericPackage';
        dependencyMatrix_ = [];
    end
    
    methods
        function obj = GenericPackage(ownerOrProcesses, outputDirectory,varargin)
            % Constructor of class Package
            
            if nargin > 0
                if(iscell(ownerOrProcesses))
                    % owner is actually a list of processes
                    obj.processes_ = ownerOrProcesses;
                    % get owner from first process
                    obj.owner_ = obj.processes_{1}.getOwner();
                else
                    % Add all processes to the generic package
                    obj.owner_ = ownerOrProcesses;
                    obj.processes_ = owner.processes_;
                end
                if(nargin > 1 && ~isempty(outputDirectory))
                    obj.outputDirectory_ = outputDirectory;
                else
                    obj.outputDirectory_ = MD.outputDirectory_;
                end
                
                nVarargin = numel(varargin);
                if nVarargin > 1 && mod(nVarargin,2)==0
                    for i=1 : 2 : nVarargin-1
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
                               
                obj.createTime_ = clock;
                
                % Make process depend on the previous process
                obj.dependencyMatrix_ = diag(ones(length(obj.processes_)-1,1),-1);
            end
        end
        function m = getDependencyMatrix(obj,i,j)   
            m = obj.dependencyMatrix_;
            if(nargin > 2)
                m = m(i,j);
            elseif(nargin > 1)
                m = m(i,:);
            end
                
        end
        function classes = getProcessClassNames(obj,index)
            classes = cellfun(@class,obj.processes_,'Unif',false);
            if(nargin > 1)
                classes = classes{index};
            end
        end
        function name = getName(obj)
            name = obj.name_;
        end
        function procConstr = getDefaultProcessConstructors(obj,index)
            procConstr = cellfun(@(x) str2func(class(x)), ...
                obj.processes_, ...
                'UniformOutput',false);
            if(nargin > 1)
                procConstr = procConstr(index);
            end
        end
    end
    
    
    methods(Static)
        function varargout = GUI(varargin)
            
            if nargin>0 && isa(varargin{1},'MovieList')
                varargout{1} = packageGUI('GenericPackage',[varargin{1}.getMovies{:}],...
                    varargin{2:end}, 'ML', varargin{1});
            else
                varargout{1} = packageGUI('GenericPackage',varargin{:});
            end

        end
        % Return the name of the package
%         function name = getName()
%             name = 'GenericPackage';
%         end
%         
    end
    
end

