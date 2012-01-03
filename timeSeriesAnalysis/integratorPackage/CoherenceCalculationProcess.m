classdef CoherenceCalculationProcess < TimeSeriesProcess
    % A concrete process for calculating coherence of sampled processes
    %
    % Sebastien Besson, Dec 2011
    
    methods (Access = public)
        
        function obj = CoherenceCalculationProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieObject'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = CoherenceCalculationProcess.getName;
                super_args{3} = @calculateMovieCoherence;
                if isempty(funParams)
                    funParams=CoherenceCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@TimeSeriesProcess(super_args{:});
        end
        
        
        function varargout = loadChannelOutput(obj,i,j,varargin)
            % Check input
            outputList={'bootstrap','avgCoh','cohCI','f'};
            ip=inputParser;
            ip.addRequired('obj');
            ip.addRequired('i',@isscalar);
            ip.addRequired('j',@isscalar);
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,i,j,varargin{:});
            output=ip.Results.output;
            if ischar(output), output={output}; end
            
            if strcmp(output{:},'bootstrap')
                s=load(obj.outFilePaths_{i,j},'avgCoh','cohCI','f');
            else
                s=load(obj.outFilePaths_{i,j},output{:});
            end
            
            for j=1:numel(output)
                if ismember(output{j},{'bootstrap'})
                    varargout{j}=s;
                else
                    varargout{j} = s.(output{j});
                end
            end
        end
        
        function output = getDrawableOutput(obj)
            output(1).name='Bootsrapped spectral density';
            output(1).var='bootstrap';
            output(1).formatData=@formatBootstrappedSpectralDensity;
            output(1).type='correlationGraph';
            output(1).defaultDisplayMethod = @(i,j)CorrelationGraphDisplay('XLabel','Frequency (Hz)',...
                'YLabel',obj.getDrawableOutputName(i,j),'Input1',obj.getInput(i).name);
            
        end
        function [label,title] = getDrawableOutputName(obj,i,j)
            if i==j, 
                label='Power spectrum (db/Hz)'; 
                title = [obj.getInput(i).name ' power spectrum'];
            else
                label='Coherence'; 
                title = [obj.getInput(i).name '/' obj.getInput(j).name ' coherence'];
            end
        end
        
    end
    
    methods (Static)
        function name =getName()
            name = 'Coherence Calculation';
        end
        function h =GUI()
            h = @coherenceCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieObject'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'coherence'];
            funParams.ProcessName=TimeSeriesProcess.getTimeSeriesProcesses;
            if isa(owner,'MovieList'),
                funParams.MovieIndex=1:numel(owner.movies_);
                winProc =cellfun(@(x) x.processes_{x.getProcessIndex('WindowingProcess',1,false)},...
                    owner.movies_,'UniformOutput',false);
                funParams.BandMin=1;
                funParams.BandMax=min(cellfun(@(x) x.nBandMax_,winProc));
                funParams.SliceIndex=cellfun(@(x) ones(x.nSliceMax_,1),winProc,'UniformOutput',false);
            else
                winProc =owner.processes_{owner.getProcessIndex('WindowingProcess',1,false)};
                funParams.BandMin=1;
                funParams.BandMax=winProc.nBandMax_;
                funParams.SliceIndex=ones(winProc.nSliceMax_,1);
            end
            funParams.nWin=8;
            funParams.window='hamming';
            funParams.noLap=.5;
            funParams.nBoot=1e4;
            funParams.alpha=.01;
        end
    end
end
% 
% function data =formatSpectralDensity(data)
% data.X=data.f;
% data.Z=10*log10(data.P);
% data=rmfield(data,{'P','f'});
% end

function data =formatBootstrappedSpectralDensity(data)
data.X=data.f;
data.Y=data.avgCoh;
data.bounds=data.cohCI;
end
