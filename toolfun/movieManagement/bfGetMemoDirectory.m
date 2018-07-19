function [ bfMemoDir ] = bfGetMemoDirectory( makeDir )
%bfGetMemoDirectory Gets the Bioformats Memoizer directory and ensures it
%exists
%
% INPUT
% makedir - boolean. If true, this will also ensure that the memo directory
%           exists. Default: true

% Mark Kittisopikul, Northwestern, 2018/01/09

    if(nargin < 1)
        makeDir = true;
    end

    bfMemoDir = [tempdir 'bioformatsMatlab'];
    
% Decided this was best dealt with upstream:
% https://github.com/openmicroscopy/bioformats/issues/3034
% id - If passed, this may modify the memo directory.
%     On windows, this will append the drive letter if given
%     if(ispc && nargin > 1 && id(2) == ':')
%         % Append Windows Drive Letter to memo directory
%         bfMemoDir = [bfMemoDir filesep id(1)];
%     end
    if(makeDir)
        w = warning;
        warning('off','MATLAB:MKDIR:DirectoryExists');
        [status] = mkdir(bfMemoDir);
        warning(w);
    end


end