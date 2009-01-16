%TRACKOBJ is a wrapper for trackCloseGapsKalman
%
% CONSTRUCTOR: obj = trackObj(coords,parameters,options)
%   IN : coords: variable containing the necessary information to construct
%                the inut argument movieInfo of trackCloseGapsKalman. Set
%                options to select the input parser.
%		parameters: (opt) Structure with fields
%                     .default : selection of default settings
%                        'standard' : standard defaults
%                        'brownian' : brownian motion only
%                     .costMatrices; gapCloseParm; kalmanFunctions
%                        optional fields to modify the default settings
%                    If parameters is not supplied, default 1 will be
%                      chosen
%		options: (opt) structure with fields
%                      .inputParser : selection of how to treat input
%                         'movieInfo' - data is already supplied as
%                            movieInfo. Parser will add standard deviations
%                            if necessary. This is the default selection if
%                            options is not supplied.
%                            inputInfo: struct with fields
%                             .fieldNames: 1-by-4 cell array of field names that
%                                should become 'xCoord', 'yCoord', 'zCoord', and
%                                'amp', respectively. Leave empty where you have no data
%                                (leaving zCoord empty will make the problem 2D)
%                                NOTE: use double-curly brackets to put a
%                                cell array into a field of a structure.
%                              .correctCentroid : 0: no correction
%                                (default), 1: correct translation, 2;
%                                correct translation and rotation
%                          'initCoord' - data is in initCoord-form, i.e.
%                             there is a field 'allCoord' with
%                             [x,y,z,sx,sy,sz], and a field amp
%                             inputInfo: struct with fields
%                               .amp - field name to use for amp. Default:
%                                 'amp'
%                               .correctCentroid -  0: no correction
%                                (default), 1: correct translation, 2;
%                                correct translation and rotation
%                      .inputInfo  : additional information for the parser
%                      .outputStyle : settings for outputStyle (see
%                          property description of outputStyle)
%                      .plotOptions : settings for plot options (see
%                          property description of plotOptions)
%                      .name : name of data set. Default ''
%                      .additionalOutput : if 1, kalmanInfoLink and errFlag
%                          are also stored, if 0, they are discarded.
%                          Default: 0;
%
%   OUT: obj: trackObject with the following properties and methods:
%
% PROPERTIES
%   movieInfo/costMatrices/gapCloseParam/kalmanFunctions/probDim
%       input for trackCloseGapsKalman
%   rawMovieInfo is the movieInfo as supplied by the user (i.e. original
%       positions)
%   tracksFinal/kalmanInfoLink/errFlag
%       output of trackCloseGapsKalman
%   results : results according to the outputStyle. results is a structure
%           with at least the fields .featureIdx and .trackInfo
%           If outputStyle is:
%  TBD!!      'separateObjects' - featureIdx is a nTracks-by-1 cell array
%                that contains a nTimepoints-by-1 list of object indices
%                with NaNs where there is no data, or a nTimepoints-by-2
%                list with [timepoint, index], if .nanList is 1 or 0,
%                respectively.
%                trackInfo is a nEvents-by-3 array with [idx,id,t], where
%                idx is the track number, id is 1,2,3, or 4 depending on
%                whether it is a birth, death, merge, or split, and t is
%                the timepoint where it happens (for merge/split t is the
%                timepoint between where the two trajectories are together
%                and where they are separate.
%                Additional options
%                   .nanList - see above
%                   .allowMerge - if 0, all tracks with merges are
%                     discarded. Default: 1
%                   .allowSplit - if 0, all tracks with splits are
%                     discarded. Default: 1
%                   .msStrategy - strategy to identify which spot is which
%                     after merge followed by split.
%   outputStyle : structure with fields
%       .name : name of the style
%       .options : additional options
%   plotOptions : structure with fields
%       timeRange,colorTime,markerType,indicateSE,axH,image,flipXY
%           image, flipXY will only be considered for 2D data
%           timeRange    : 2-element row vector indicating time range
%                           to plot.  Optional. Default: whole movie.
%           colorTime    : String with the following options:
%                           -'1' if time is to be color-coded (green in the
%                           beginning, blue in the middle, red in the end).
%                           -'k', 'b', 'r', etc. if all tracks are in black,
%                           blue, red, etc.
%                           Optional. Default: 'k'.
%           markerType   : String indicating marker type for plotting.
%                           Only used if colorTime is not '1'.
%                           Optional. Default: 'none'.
%           indicateSE   : 1 if track starts and ends are to be indicated
%                           with circles and squares, respectively; 0
%                           otherwise. Optional. Default: 1.
%           axH          : if non-empty, data will be plotted into axes
%                           specified by axH. If 1, a new figure will be
%                           opened. If 0, plot will be into current axes.
%                           Default: 0
%           image        : An image that the tracks will be overlaid on if
%                           newFigure=1. It will be ignored if newFigure=0.
%                           Optional. Default: no image
%           flipXY       : 1 if x and y coord should be flipped for
%                           plotting. Optional. Default: 0.
%           useRaw       : 1 if raw (unaligned) coordinates should be used
%                           for plotting. Optional. Default: 0;
%           beforeAfter  : 2-element vector with # of timepoints before and
%                           after current timepoint that is to be plotted.
%                           Used with method plotFrame
%
%    name : name of data set
%
% METHODS
%
%   trackObj = trackObj(coords,parameters,options)
%       Constructor (see above)
%
%   [tracksFinal,kalmanInfoLink,errFlag] = run(trackObj,verbose,saveResults)
%       This is the main tracking function
%       verbose      : 1 to show calculation progress, 0 otherwise.
%                      Optional. Default: 1.
%       saveResults  : 0 if no saving is requested. Default
%                      If saving is requested, structure with fields:
%           .dir          : Directory where results should be saved.
%                           Optional. Default: current directory.
%           .filename     : Name of file where results should be saved.
%                      Or []. Default: trackedFeatures in directory
%                      where run is initiated.
%                      Whole structure optional.
%       For description of output, please see tracksCloseGapsKalman
%
%    plot(trackObj,plotOpt)
%       Plots the track result. Uses options defined in obj.plotOptions,
%       plotOpt is a structure of the same form as plotOptions and allows
%       to overrule plotOptions.
%
%    plotFrame(trackObj,t,ah)
%       plots the track around frame t into the axes specified by ah
%
% REMARKS trackObj is a handle class. Thus, it will be passed by reference!
%         Note this is still work in progress.
%
% created with MATLAB ver.: 7.7.0.2162 (R2008b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 14-Aug-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef trackObj<handle
    properties
        % input
        movieInfo
        costMatrices
        gapCloseParam
        kalmanFunctions
        probDim
        % output
        tracksFinal = [];
        kalmanInfoLink = [];
        errFlag = [];
        % additional properties
        outputStyle = struct('name','separateObjects','options','');
        plotOptions = struct('timeRange',[],'colorTime','1',...
            'markerType','none','indicateSE',1,'axH',0,...
            'image',[],'flipXY',0,'useRaw',0);
        name = '';
        additionalOutput  = false;
        nTimepoints = 0;
        rawMovieInfo % raw coordinates (before alignment)

    end % static props
    properties (Dependent)
        results % property depending on outputStyle
    end % dependent properties
    properties (Hidden)
        % set here default options for outputStyles. For every name, there
        % need to be defaults for all corresponding options
        outputStyleDefOpt = struct('separateObjects',...
            struct('nanList',true,'allowMerge',true,'allowSplit',true,...
            'msStrategy','ampDist')...
            );
        % defaults for inputInfo
        % correctCentroid: 0 - none; 1- only centroid; 2 - centroid and rot
        % findNames: check for fieldnames
        inputInfoDefOpt = struct('movieInfo',...
            struct('fieldNames',{{'xCoord','yCoord','zCoord','amp'}},...
            'correctCentroid',0,'findNames',0),...
            'initCoord',...
            struct('correctCentroid',0));
        
    end % hidden properties
    methods
        %================
        %% CONSTRUCTOR
        %================
        function obj = trackObj(coords,parameters,options)
            % allow empty object, in case we want to run a batch of tracks
            if nargin < 1 || isempty(coords)
                obj.movieInfo = [];
                return
            end
            
            % set defaults
            if nargin < 2 || isempty(parameters) || ~isfield(parameters,'default')
                parameters.default = 'brownian';
            end
            if nargin < 3 || isempty(options) || ~isfield(options,'inputParser')
                options.inputParser = 'movieInfo';
            end
            if ~isfield(options,'inputInfo')
                options.inputInfo = [];
            end
            
            % set options
            options = obj.setOptions(options);
            
            % parse input to set movieInfo. Return movieInfo explicitly so
            % that we can call parseInputs from itself, if necessary
            obj.movieInfo = obj.parseInputs(coords,options);
            
            % set parameters
            obj.setParameters(parameters);
            
            
            
        end
        %===================
        %% RUN
        %===================
        function [tracksFinal,kalmanInfoLink,errFlag] = run(obj,verbose,saveResults)
            
            %---- test input
            nObj = numel(obj);
            if nObj > 1 && nargout > 0
                error('cannot return output arguments with multiple objects')
            end
            if nargin < 2 || isempty(verbose);
                verbose = 1;
            end
            if nargin < 3 || isempty(saveResults)
                saveResults = zeros(nObj,1);
            else
                nSaveResults = length(saveResults);
                if nSaveResults == 1 && nObj>1
                    saveResults = repmat(saveResults,nObj,1);
                end
                if nSaveResults ~= nObj
                    error('for multiple objects, multiple save options must be provided')
                end
            end
            
            %---- track
            for iObj = 1:nObj
                [obj(iObj).tracksFinal,obj(iObj).kalmanInfoLink,obj(iObj).errFlag]...
                    = trackCloseGapsKalman(obj(iObj).movieInfo,...
                    obj(iObj).costMatrices,obj(iObj).gapCloseParam,...
                    obj(iObj).kalmanFunctions,obj(iObj).probDim,saveResults(iObj),verbose);
            end
            
            %---- return results if requested
            if nargout > 0
                % nObj equals to 1
                tracksFinal = obj.tracksFinal;
                kalmanInfoLink = obj.kalmanInfoLink;
                errFlag = obj.errFlag;
            end
            % discard kalmanInfoLink if unnecessary
            for iObj = 1:nObj
                if obj(iObj).additionalOutput == 0
                    obj(iObj).kalmanInfoLink = [];
                    obj(iObj).errFlag = [];
                end
            end
            
        end % run
        %======================
        %% PLOT
        %======================
        function plot(obj,plotOpt)
            % todo: allow direct input to overwrite plotOptions
            % todo: set callback (don't forget to set userData in
            % plotTracks)
            % loop through objects and plot
            
            
            
            nObj = numel(obj);
            for iObj = 1:nObj
                
                % read plotOptions individually for each object
                if nargin < 2 || isempty(plotOpt)
                    plotOptions = obj(iObj).plotOptions;
                else
                    plotOptions = obj(iObj).plotOptions;
                    % overwrite plotOptions
                    for fn = fieldnames(plotOpt)'
                        plotOptions.(fn{1}) = plotOpt.(fn{1});
                    end
                end
                
                
                tmpAxH = plotOptions.axH;
                if tmpAxH == 0 || (nObj>1 && tmpAxH == 1);
                    % create figure using the name that has been provided
                    if ~isempty(obj(iObj).name)
                        figure('Name',obj(iObj).name);
                    else
                        figure('Name',sprintf('Data set %i',iObj));
                    end
                    plotOptions.plotOptions.axH = 0;
                    axes(gca);
                elseif ishandle(tmpAxH) && strcmp(get(tmpAxH,'type'),'axes')
                    plotOptions.axH = 0;
                    axes(tmpAxH);
                end
                switch obj(iObj).probDim
                    case 2
                        plotTracks2D(obj(iObj).tracksFinal,...
                            plotOptions.timeRange,...
                            plotOptions.colorTime,...
                            plotOptions.markerType,...
                            plotOptions.indicateSE,...
                            plotOptions.axH,...
                            plotOptions.image,...
                            plotOptions.flipXY,0);
                    case 3
                        plotTracks3D(obj(iObj).tracksFinal,...
                            plotOptions.timeRange,...
                            plotOptions.colorTime,...
                            plotOptions.markerType,...
                            plotOptions.indicateSE,...
                            plotOptions.axH);
                end
                % update plotOptions with axH
                obj(iObj).plotOptions.axH = tmpAxH;
            end
        end % plot
        %==========================
        %% PLOTFRAME
        %==========================
        function plotFrame(obj,t,axH)
            if nargin < 3 || isempty(t) || isempty(axH)
                error('please supply time and axes handle for plotting the current frame')
            end
            
        end
        %==========================
        %% GET.RESULTS
        %==========================
        % if there is a list of objects, return multiple outputs, as we
        % would expect of, e.g. a structure
        function varargout = get.results(obj)
            % count objects
            nObj = numel(obj);
            % assign default output
            [varargout{1:nObj}] = deal(struct('featureIdx',[],...
                'trackInfo',[]));
            
            % loop objects.
            for iObj = 1:numel(obj)
                %Since obj is passed by reference, we can copy
                % to make our lives easier
                cObj = obj(iObj);
                % only read if there is anything to do
                if ~isempty(cObj.tracksFinal)
                    % switch according to outputStyle
                    switch cObj.outputStyle.name
                        case 'separateObjects'
                            % merges have to be resolved into multiple
                            % tracks.
                            
                            % pass once through the list to get the number
                            % of segments and to remove merges/splits/ms if
                            % necessary
                            %                             nTracks = 0;
                            %                             trackList = cObj.tracksFinal;
                            %                             for i=length(trackList):-1:1
                            %                                 % check for merge/split
                            %                                 nSegments = size(trackList(i).tracksFeatIndxCG,1);
                            %                                 if nSegments > 1 && ...
                            %                                         any(tracksFinal(iTrack).seqOfEvents(:,2)==2 & ...
                            %                                         ~isnan(tracksFinal(iTrack).seqOfEvents(:,4)))
                            %                                     % bad track
                            %                                     nSegments = 0;
                            %                                 end
                            %                             end
                            
                        otherwise
                            error('unrecognized output style %s',cObj.outputStyle.name)
                    end
                end
            end
        
        end
        %% Get Method rawMovieInfo
        function out = get.rawMovieInfo(obj)
            % return movieInfo unless correctCoords has written
            % rawMovieInfo
            if isempty(obj.rawMovieInfo)
                out = obj.movieinfo;
            else
                out = obj.rawMovieInfo;
            end
        end
    end % public methods
    methods (Hidden)
        %===========================
        %% PARSE INPUTS
        %===========================
        function movieInfo = parseInputs(obj,coord,options)
            switch options.inputParser
                case 'movieInfo'
                    % set default field names
                    movieInfoFields = options.inputInfo.fieldNames;
                    if ~options.inputInfo.findNames
                        fieldList = options.inputInfo(:);
                    else
                        % if no explicit list: find fieldnames from coord
                        fieldList = fieldnames(coord);
                        nameIdx = ismember(movieInfoFields',fieldList);
                        if nameIdx(1)==0
                            error('you need at least to supply x coords')
                        end
                        fieldList = movieInfoFields;
                        [fieldList{~nameIdx}] = deal('');
                    end
                    if length(fieldList) == 3 || isempty(fieldList{3})
                        % 2D problem. remove zCoord
                        obj.probDim = 2;
                        movieInfoFields(3) = [];
                        if isempty(fieldList{3})
                            fieldList(3) = [];
                        end
                    else
                        obj.probDim = 3;
                    end
                    % create movieInfo
                    nTimepoints = length(coord);
                    tmp = movieInfoFields;
                    tmp{2,end} = [];
                    movieInfo(1:nTimepoints) = struct(tmp{:});
                    
                    % loop to fill movieInfo
                    for iName = 1:1+obj.probDim
                        for t = 1:nTimepoints
                            if ~isempty(fieldList(iName))
                                movieInfo(t).(movieInfoFields{iName}) = ...
                                    coords(t).(fieldList{iName});
                                if size(movieInfo(t).(movieInfoFields{iName}),2) == 1
                                    movieInfo(t).(movieInfoFields{iName})(:,2) = 1;
                                end
                            else
                                % set the un-supplied data to 1. This makes
                                % only really sense for amp, though
                                movieInfo(t).(movieInfoFields{iName}) = ...
                                    ones(size(movieInfo(t).xCoord));
                            end
                        end
                    end
                    
                    % correct centroid
                    correctCentroid;
                    
                case 'initCoord'
                    % check options
                    if ~isempty(options.inputInfo) && isfield(options.inputInfo,'amp')
                        ampName = options.inputInfo.amp;
                    else
                        ampName = 'amp';
                    end
                    
                    % problem dimension is automatically 3
                    obj.probDim = 3;
                    % create movieInfo
                    nTimepoints = length(coord);
                    movieInfo(1:nTimepoints) = struct('xCoord',[],...
                        'yCoord',[],'zCoord',[],'amp',[]);
                    
                    goodTimes = find(arrayfun(@(x)(~isempty(x.allCoord)),...
                        coord));
                    
                    if ~isempty(goodTimes)
                        % only do the rest if nonempty input
                        
                        for t = goodTimes'
                            movieInfo(t).xCoord = coord(t).allCoord(:,[1 4]);
                            movieInfo(t).yCoord = coord(t).allCoord(:,[2 5]);
                            movieInfo(t).zCoord = coord(t).allCoord(:,[3 6]);
                            movieInfo(t).amp = coord(t).(ampName);
                        end
                        
                        correctCentroid;
                        
                    end
                    
            end % switch inputParser
            
            obj.nTimepoints = length(movieInfo);
            
            %--- nested function correctCentroid
            function correctCentroid
                % correct Centroid
                        if options.inputInfo.correctCentroid > 0
                            % store rawMovieInfo
                            obj.rawMovieInfo = movieInfo;
                            
                            % if not correct rotation, co is -1, otherwise inf
                            if options.inputInfo.correctCentroid == 1;
                                co = -1;
                            else
                                co = Inf;
                            end
                            [movieInfo,obj.outputStyle.centroids,...
                                obj.outputStyle.rotMats] = alignCoords(...
                                movieInfo,{'xCoord','yCoord','zCoord','amp'},...
                                co,obj.probDim);
                        end
            end
        end % parseInputs
        %===========================
        %% SET PARAMETERS
        %===========================
        function setParameters(obj,parameters)
            % read defaults
            setDefaultParameters(obj,parameters.default);
            
            % go through parameters and overwrite defaults. Blindly step
            % through fieldnames and overwrite whatever was in the input.
            % This requires that the user made sensible choices.
            for pn = fieldnames(parameters)'
                if ~strcmp(pn{1},'default')
                    for k = 1:length(parameters.(pn{1}))
                        for ppn = fieldnames(parameters.(pn{1})(k))
                            if isstruct(parameters.(pn{1})(k).(ppn{1}))
                                for pppn = fieldnames(parameters.(pn{1})(k).(ppn{1}))'
                                    try
                                        obj.(pn{1})(k).(ppn{1}).(pppn{1}) = ...
                                            parameters.(pn{1})(k).(ppn{1}).(pppn{1});
                                    catch
                                        error(...
                                            'parameters.%s(%i).%s.%s is not a recognized option',...
                                            pn{1},k,ppn{1},pppn{1})
                                    end
                                end
                            else
                                try
                                    obj.(pn{1})(k).(ppn{1}) = ...
                                        parameters.(pn{1})(k).(ppn{1});
                                catch
                                    error(...
                                        'parameters.%s(%i).%s is not a recognized option',...
                                        pn{1},k,ppn{1})
                                end
                            end
                        end
                    end
                end
            end
            
            
        end % set parameters
        %===========================
        %% SET OPTIONS
        %===========================
        function options = setOptions(obj,options)
            %		options: (opt) structure with fields
            %                      .inputParser : selection of how to treat input
            %                         'movieInfo' - data is already supplied as
            %                            movieInfo. Parser will add standard deviations
            %                            if necessary. This is the default selection if
            %                            options is not supplied.
            %                      .inputInfo  : additional information for the parser
            %                      .outputStyle : settings for outputStyle (see
            %                          property description of outputStyle)
            %                      .plotOptions : settings for plot options (see
            %                          property description of plotOptions)
            %                      .name : name of data set. Default ''
            %                      .additionalOutput : if 1, kalmanInfoLink and errFlag
            %                          are also stored, if 0, they are discarded.
            %                          Default: 0;
            % options will never be empty (there is at least the default
            % input parser. Therefore, just test for fieldnames
            
            % outputStyle
            if isfield(options,'outputStyle')
                % set defaults
                if isfield(options.outputStyle,'name')
                    obj.outputStyle.name = options.outputStyle.name;
                end
                obj.outputStyle.options = obj.outputStyleDefOpt.(...
                    obj.outputStyle.name);
                % overwrite default with additional options
                for fn = fieldnames(options.outputStyle)'
                    if strcmp(fn{1},'options')
                        for ofn = fieldnames(options.outputStyle.options)'
                            obj.outputStyle.options.(ofn{1}) = ...
                                options.outputStyle.options.(ofn{1});
                        end
                    else
                        obj.outputStyle.(fn{1}) = options.outputStyle.(fn{1});
                    end
                end
            end
            
            % inputInfo
            if isempty(options.inputInfo)
                options.inputInfo = obj.inputInfoDefOpt.(options.inputParser);
            else
                tmp = options.inputInfo;
                options.inputInfo = obj.inputInfoDefOpt.(options.inputParser);
                for fn = fieldnames(tmp)'
                    options.inputInfo.(fn{1}) = tmp.(fn{1});
                end
                if strcmp(options.inputParser,'movieInfo') && ...
                        (~isfield(tmp,'fieldNames') || isempty(tmp,'fieldNames'))
                    options.inputInfo.findNames = true;
                end
            end
            
            if isfield(options,'plotOptions')
                for fn = fieldnames(options.plotOptions)'
                    obj.plotOptions.(fn{1}) = options.plotOptions.(fn{1});
                end
            end
            if isfield(options,'name')
                obj.name = options.name;
            end
            if isfield(options,'additionalOutput')
                obj.additionalOutput = options.additionalOutput;
            end
        end % set options
        %===========================
        %% DEFAULTPARAMETERS
        %===========================
        function setDefaultParameters(obj,default)
            switch lower(default)
                case 'brownian'
                    % general gap closing parameters
                    obj.gapCloseParam.timeWindow = 5; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
                    obj.gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 0 otherwise.
                    obj.gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.
                    
                    % cost matrix for frame-to-frame linking
                    
                    %function name
                    obj.costMatrices(1).funcName = 'costMatLinearMotionLink';
                    
                    %parameters
                    
                    parameters.linearMotion = 0; %no linear motion
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
                    parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
                    
                    parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
                    
                    obj.costMatrices(1).parameters = parameters;
                    clear parameters
                    
                    % cost matrix for gap closing
                    
                    %function name
                    obj.costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';
                    
                    %parameters
                    
                    %needed all the time
                    parameters.linearMotion = 0;%1; %use linear motion Kalman filter.
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius.
                    parameters.brownStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
                    parameters.timeReachConfB = 2; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).
                    
                    parameters.ampRatioLimit = [0.5 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
                    
                    parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
                    
                    parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
                    
                    parameters.linStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
                    parameters.timeReachConfL = obj.gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
                    parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
                    
                    obj.costMatrices(2).parameters = parameters;
                    clear parameters
                    
                    % Kalman filter function names
                    
                    obj.kalmanFunctions.reserveMem = 'kalmanResMemLM';
                    obj.kalmanFunctions.initialize = 'kalmanInitLinearMotion';
                    obj.kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
                    obj.kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
                case 'standard'
                    % general gap closing parameters
                    obj.gapCloseParam.timeWindow = 5; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
                    obj.gapCloseParam.mergeSplit = 1; %1 if merging and splitting are to be considered, 0 otherwise.
                    obj.gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.
                    
                    % cost matrix for frame-to-frame linking
                    
                    %function name
                    obj.costMatrices(1).funcName = 'costMatLinearMotionLink';
                    
                    %parameters
                    
                    parameters.linearMotion = 1; %use linear motion Kalman filter.
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
                    parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
                    
                    parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
                    
                    obj.costMatrices(1).parameters = parameters;
                    clear parameters
                    
                    % cost matrix for gap closing
                    
                    %function name
                    obj.costMatrices(2).funcName = 'costMatLinearMotionCloseGaps';
                    
                    %parameters
                    
                    %needed all the time
                    parameters.linearMotion = 1;%1; %use linear motion Kalman filter.
                    
                    parameters.minSearchRadius = 1.5;%1.5; %minimum allowed search radius.
                    parameters.maxSearchRadius = 20.0;%1.5; %maximum allowed search radius.
                    parameters.brownStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
                    parameters.timeReachConfB = 2; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).
                    
                    parameters.ampRatioLimit = [0.5 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
                    
                    parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
                    
                    parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
                    parameters.nnWindow = obj.gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
                    
                    parameters.linStdMult = 3*ones(obj.gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
                    parameters.timeReachConfL = obj.gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
                    parameters.maxAngleVV = 45; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
                    
                    obj.costMatrices(2).parameters = parameters;
                    clear parameters
                    
                    % Kalman filter function names
                    
                    obj.kalmanFunctions.reserveMem = 'kalmanResMemLM';
                    obj.kalmanFunctions.initialize = 'kalmanInitLinearMotion';
                    obj.kalmanFunctions.calcGain = 'kalmanGainLinearMotion';
                    obj.kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';
                otherwise
                    error('%s not implemented yet',default)
            end
        end
       
    end % hidden methods
    
end