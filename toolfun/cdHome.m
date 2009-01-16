function path = cdHome(opt)
%CDHOME changes directories to HOME and/or returns the HOME path
%
% SYNOPSIS: path = cdHome(opt)
%
% INPUT opt: optional options
%		0: only return path
%		1: return path and change directory
%       2: return old path and change directory
%
% OUTPUT path: MATLABHOME or HOME
%
% REMARKS matlabhome/home is $HOME/matlab
%
% created with MATLAB ver.: 7.7.0.471 (R2008b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 31-Dec-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TEST INPUT
if nargin < 1 || isempty(opt)
    opt = 1;
end


%% FIND HOME
switch opt
    case {0,1}
        path=getenv('MATLABHOME');
        if isempty(path)
            path=getenv('HOME');
            if isempty(path)
                error('Neither MATLABHOME nor HOME are defined')
            end
        end
        % append matlab
        try
            mPath = fullfile(path,'matlab');
            op = cd(mPath);
            cd(op);
            path = mPath;
        end
            
    case 2
        path = pwd;
        cdHome(1);
    otherwise
        error('opt %i not implemented yet',opt)
end

%% CHANGE DIRECTORY
switch opt
    case {0,2}
        % do nothing
    case 1
        % change directory
        cd(path)
end
