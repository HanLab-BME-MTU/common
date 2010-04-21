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
        processes_ % object array to store process objects
        packages_ % object array to store package objects
        movieDataPath_
        movieDataFileName_
        notes_
    end

    methods
        function obj = MovieData(channelPath, pixelSize, timeInterval,...
                movieDataPath, movieDataFileName, notes)
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
        %% MovieData Sanity Check
        function isValid = sanityCheck(obj, movieDataPath, ...
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
            % Check saved user input is numeric 
            if ~isnumeric(obj.pixelSize_ ) || ...
                                     ~isnumeric(obj.timeInterval_)
               % Exception: 
               error('LCCB:SanMD:WrongMInfo',...
                   'For unknown reason, the value of pixel size or time interval is invalid\n\n');
            end
            % I. Check if the path and filename stored in the movieData are the same
            % than the ones provided in argument. They can differ if the movieData
            % file has been rename, move or copy to another location.
            if  ~strcmp(obj.movieDataPath_, movieDataPath)
               obj.movieDataPath_ = movieDataPath;
            end
            if  ~strcmp(obj.movieDataFileName_, movieDataFileName)
               obj.movieDataFileName_ = movieDataFileName; 
            end
            
            % Dig into every channel path
            imSize = zeros(2, length(obj.channelPath_));
            nFrames = zeros(length(obj.channelPath_),1);
            fExt = {'tif','TIF','STK','bmp','BMP','jpg','m'};
            for i=1:length(obj.channelPath_)
                % II. Check that every channelPath exists
                if isempty( dir(obj.channelPath_{i}) )
                    % Exception: 
                    error('LCCB:SanMD:NoPath',...
                        ['One or more channel paths do not exist. '...
                        'Please make sure the channel path/paths are corrent\n\n']);
                end
                % Check the number of file extensions
                fileNames = [ ];
                nofExt = 0;
                for j = 1: length(fExt)
                    tempfileNames = dir([obj.channelPath_{i}, ...
                                                 filesep, '*.', fExt{j}]);
                    if ~isempty(tempfileNames)
                        nofExt = nofExt + 1;
                    end
                    fileNames = vertcat(fileNames, tempfileNames);                  %#ok<AGROW>
                end
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
            
            isValid =true;

        end
        
        %%
        function addChannelPath(obj, newpath)
            % Add a new channel path to the object
            if isempty(obj.channelPath_)
               obj.channelPath_{1} = newpath; 
            else
                obj.channelPath_{end+1} = newpath;
            end 
        end
        %% Functions to manipulate process object array
        function addProcesses(obj, newprocess)
            % Add a process to the processes_ array
            if isempty(obj.processes_)
                obj.processes_ = newprocess;
            else
                obj.processes_(end+1) = newprocess;
            end
        end
        function val = getProcess(obj, i)
            % Get the i th process in processes_ array
            if i > length(obj.processes_)
                error('Index of process exceeds dimension of processes_ array');
            else
                val = obj.processes_(i);
            end
        end
        function setProcess(obj, i, process)
            % Set the i th process in processes_ array
            if i > length(obj.processes_)
                error('Index of process exceeds dimension of processes_ array');
            else
                delete(obj.processes_(i));
                obj.processes_(i) = process;
            end
        end
        %% Functions to manipulate package object array
        function addPackages(obj, newpackage)
            % Add a package to the packages_ array
            if isempty(obj.packages_)
                obj.packages_ = newpackage;
            else
                obj.packages_(end+1) = newpackage;
            end
        end
        function val = getPackage(obj, i)
            % Get the i th package in packages_ array
            if i > length(obj.packages_)
                error('Index of package exceeds dimension of packages_ array');
            else
                val = obj.packages_(i);
            end
        end
        function setPackage(obj, i, package)
            % Set the i th package in packages_ array
            if i > length(obj.packages_)
                error('Index of package exceeds dimension of packages_ array');
            else
                delete(obj.packages_(i));
                obj.packages_(i) = package;
            end
        end

    end

    
end