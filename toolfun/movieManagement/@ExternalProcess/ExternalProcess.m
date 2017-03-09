classdef ExternalProcess < NonSingularProcess
    % ExternalProcess A minimalist implementation of a Process that can run
    % an arbitrary function with the ExternalProcess handle as an argument
    %
    % Constructor
    % -----------
    % ExternalProcess(owner,name,fun)
    % owner: A MovieObject such as MovieData
    % name:  (Optional) A char string describing the process
    % fun:   (Optional) A function handle taking the process as an
    %                   argument
    %
    % Static Methods
    % .getName()          : 'Dummy process'
    % .getDefaultParams() : struct()
    %
    % Example
    % -------
    % MD = MovieData;
    % process = ExternalProcess(MD,'Say something',@(p) disp(p.getParameters().text));
    % process.setParameters(struct('text','Hello world!'))
    % process2 = ExternalProcess(MD,'Say something',@(p) disp(p.getParameters().text));
    % process2.setParameters(struct('text','Good bye!'))
    % MD.addProcess(process)
    % MD.addProcess(process2);
    % cellfun(@run,MD.processes_);
    
    methods(Access = public)
        
        function obj = ExternalProcess(owner, varargin)
            %ExternalProcess(owner,name,fun)
            %owner: A MovieObject such as MovieData
            %name:  (Optional) A char string describing the process
            %fun:   (Optional) A function handle taking the process as an
            %                  argument
            
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('name',ExternalProcess.getName(),@ischar);
            ip.addOptional('fun',@(x) x,@(f) validateattributes(f,{'function_handle','char'},{}));
            ip.parse(owner,varargin{:});
            
            % Constructor of the ExternalProcess
            super_args{1} = owner;
            super_args{2} = ip.Results.name;
            obj = obj@NonSingularProcess(super_args{:});
            obj.funName_ = ip.Results.fun;
            obj.funParams_ = obj.getDefaultParams(owner,varargin{:});
            
        end
    end
    methods (Static)
        function name = getName()
            name = 'Dummy process';
        end
        function funParams = getDefaultParams(varargin)
            funParams = struct();
        end
        function func = GUI(varargin)
            func = @cliGUI;
        end
        
    end
end
