classdef  MovieData < handle
% This is the movie management class definition
    properties (SetAccess = private, GetAccess = public)
    % SetAccess = private - cannot change the values of variables outside object
    % GetAccess = public - can get the values of variables outside object without
    % definging accessor functions
        nFrames_ 
        imSize_
        channelPath_
        pixelSize_
        timeInterval_
        processes_ = {}; % Cell array to store process objects
        packages_ = {};% Cell array to store package objects
        movieDataPath_
        movieDataFileName_
        notes_
        
        crtPackage_ = [ ]; % Handle of current package
        crtProcess_ = [ ]; % Handle of current process
    end

    methods
        function obj = MovieData(channelPath, pixelSize, timeInterval,...
                movieDataPath, movieDataFileName, notes) % throws exception
            % MovieManager construntor
            if nargin > 0
                
                % Input validation: channelPath
                if isempty(channelPath)
                   % Exception
                   error('LCCB:MovieData:NoPath',...
                       'Please provide at least one channel path\n\n');
                else
                    obj.channelPath_ = channelPath;
                end
                % Input validation pixelSize and timeInterval
                if isempty(pixelSize) || isempty(timeInterval)
                   % Exception
                   error('LCCB:MovieData:NoInput',...
                       'Please provide pixel size and time interval values\n\n');
                elseif isnumeric(pixelSize) && isnumeric(timeInterval)
                    obj.pixelSize_ = pixelSize;
                    obj.timeInterval_ = timeInterval;
                elseif isnan(str2double(pixelSize)) || ...
                                     isnan(str2double(timeInterval))
                        % Exception
                    error('LCCB:MovieData:WrongMInfo',...
                        'Movie information is invalid\n\n');
                elseif str2double(pixelSize)<=0 || ...
                        str2double(timeInterval)<=0 
                    error('LCCB:MovieData:NegaMInfo',...
                        'Movie information cannot be negative\n\n');
                else
                    obj.pixelSize_ = str2double(pixelSize);
                    obj.timeInterval_ = str2double(timeInterval);
                end
                obj.movieDataPath_ = movieDataPath;
                obj.movieDataFileName_ = movieDataFileName;
                obj.notes_ = notes;
                
            end
        end
        % MovieData Sanity Check
        function sanityCheck(obj, movieDataPath, ...
                                                movieDataFileName, full)
        % 1. Sanity check (user input, input channels, image files)
        % 2. Assignments to 4 properties:
        %       movieDataPath_
        %       movieDataFileName_
        %       nFrames_
        %       imSize_
            
            if nargin < 4 || isempty(full)
               full = false; 
            end
            % I. Check if the path and filename stored in the movieData are the same
            % as the ones provided in argument. They can differ if the movieData
            % MAT file has been rename, move or copy to another location.
            if  ~strcmp(obj.movieDataPath_, movieDataPath)
               obj.movieDataPath_ = movieDataPath;
            end
            if  ~strcmp(obj.movieDataFileName_, movieDataFileName)
               obj.movieDataFileName_ = movieDataFileName; 
            end
            
            % Dig into every channel path
            imSize = zeros(2, length(obj.channelPath_));
            nFrames = zeros(length(obj.channelPath_),1);
            for i=1:length(obj.channelPath_)
                % II. Check that every channelPath exists
                if isempty( dir(obj.channelPath_{i}) )
                    % Exception: 
                    error('LCCB:SanMD:NoPath',...
                        ['One or more channel paths do not exist. '...
                        'Please make sure the channel path/paths are corrent\n\n']);
                end
                % Check the number of file extensions
                [fileNames nofExt] = imDir(obj.channelPath_{i},true);
                switch nofExt
                    case 0
                        % ToDo: Exception: No proper image files are detected
                        error('LCCB:SanMD:NoImFile',...
                            ['No proper image files are detected in:\n'...
                                                    obj.channelPath_{i}]);
                    case 1
                        % Number of image files in each channel
                        nFrames(i) = length(fileNames); 
                    otherwise
                        % Exception: More than one type of image 
                        % files are in the current specific channel
                        error('LCCB:SanMD:ImgEx',...
                            ['More than one type of image files are found in:\n' ...
                              obj.channelPath_{i},...
                              '\n\nPlease make sure all images have the same size.\n\n']);
                end
                % III. Check the consistency of image size within single
                % channel
                imSize2 = zeros(2, length(fileNames));
                for j = 1: length(fileNames)
                    imInfo = imfinfo([ obj.channelPath_{i}, ...
                                             filesep, fileNames(j).name]);
                    imSize2(1,j) = imInfo.Width;
                    imSize2(2,j) = imInfo.Height;
                end
                if max(imSize2(1,:))~=min(imSize2(1,:)) || ...
                        max(imSize2(2,:))~=min(imSize2(2,:))
                    % Exception: Image sizes are inconsistent in the
                    % current channel. 
                    error('LCCB:SanMD:ImgSizeIn',...
                        ['Image sizes are inconsistent in:\n',... 
                         obj.channelPath_{i},...
                         '\n\nPlease make sure all the images have the same size.\n\n']);
                end
                imSize(1,i) = imSize2(1);
                imSize(2,i) = imSize2(2);
            end
            % VI. Check that each channelPath contains the same number of image files.
            if max(nFrames) ~= min(nFrames)
                % Exception: 
                error('LCCB:SanMD:ImgNum',...
                    'Different number of frames are detected in channels\n\n');
            end
            % V. Check the consistency of image size between different
            % channels
            if max(imSize(1,:))~=min(imSize(1,:)) ||...
                    max(imSize(2,:))~=min(imSize(2,:))
               % Exception: Image sizes are inconsistent between 
               % different channels 
               error('LCCB:SanMD:ImgSizeOut',...
                   'Image sizes are inconsistent in different channels\n\n');
            end
		
		% VI. Assign imSize_ and nFrames_;     
            if ~isempty(obj.nFrames_)
                if (obj.nFrames_ ~= nFrames(1))
                % Exception: the number of frames stroed in this movieData
                % is different from the number of images in channel
                % directories.
                    error('LCCB:SanMD:ReloadImgNum',...
                        'Record shows the number of frames has changed in this movie\n\n');
                end
            else
                obj.nFrames_ = nFrames(1);
            end
            if ~isempty(obj.imSize_)
                if (obj.imSize_(1) ~= imSize(1) || obj.imSize_(2) ~=imSize(2))
                % Exception: the image size stored in this movieData is
                % different from the image size of the image files
                    error('LCCB:SanMD:ReloadImgSize',...
                        'Record shows image size has changed in this movie\n\n');
                end
            else
                obj.imSize_ = imSize(:,1);
            end

            if full
            % Call the sanity check of every process and package
            end

        end
        
        %
        function addChannelPath(obj, newpath)
            % Add a new channel path to the object
            obj.channelPath_ = horzcat(obj.channelPath_, {newpath});
        end
        % Functions to manipulate process object array
        function addProcess(obj, newprocess)
            % Add a process to the processes_ array
            obj.processes_ = horzcat(obj.processes_, {newprocess});
        end
        function setProcess(obj, i, process)
            % Set the i th process in processes_ array
            assert( i > length(obj.processes_),...
                'Index of process exceeds dimension of processes_ array');
                delete(obj.processes_{i});
                obj.processes_{i} = process;
        end
        function setCrtProcess(obj, process)
           % Set the current process to 'process'
           obj.crtProcess_ = process;
        end
        % Functions to manipulate package object array
        function addPackage(obj, newpackage)
            % Add a package to the packages_ array
            obj.packages_ = horzcat(obj.packages_ , {newpackage});
        end

        function setPackage(obj, i, package)
            % Set the i th package in packages_ array
            assert( i > length(obj.packages_), ...
                'Index of package exceeds dimension of packages_ array');
                delete(obj.packages_{i});
                obj.packages_{i} = package;
        end
        function setCrtPackage(obj, package)
           % Set the current package to 'package'
           obj.crtPackage_ = package;
        end

    end

    
end