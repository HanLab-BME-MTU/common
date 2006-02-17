function gencodetemplate(fxnname)
% GENCODETEMPLATE creates a template for new functions
%
% DESCRIPTION: This function generates a template based on user-defined
% variables for new functions.  The point is to improve the documentation
% process during the software development phase.
%
% SYNOPSIS  gencodetemplate(fxnname)
%
% INPUT  fxnname (optional): name of the function that should be created
%
% OUTPUT The output is an m-file function filled in with user-defined
%         variables.
%
% MATLAB VERSION (originally written on): 7.1.0.246 (R14) Service Pack 3
%
% USERNAME: jkunken
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===============
% CHECK INPUT
%===============

% if no input argument, assign empty to functionname
if nargin < 1 || isempty(fxnname) || ~ischar(fxnname)
    fxnname = [];
end

%==============


%================
% GATHER DATA
%================

% data-gathering stage begins below
% getenv() and version are used to derive environmental variables

% find username. This will return empty on Linux, but we have a line for
% the username in the input dialogue, anyway.
username = getenv('username');

% get version, OS, date
vers    = version;
os      = getenv('OS');
datetoday = date;

% ask for the rest of the input with inputdlg

% set up inputdlg
inPrompt = {'Username',...
    'Function name',...
    'Description: ''FUNCTIONNAME does ...'' (captitalized function name)',...
    'Synopsis: ''[output1, output2] = functionname(input1, input2)',...
    'Description of input arguments (use \n for line breaks)',...
    'Description of output arguments (use \n for line breaks)'};
inTitle = 'Please describe your function';
numLines = repmat([1,100],6,1);
% assign defaultAnswer
[defaultAnswer{1:6}] = deal('');
% assign username in default answers
defaultAnswer{1} = username;
% if we have a function name already, all will be simpler
if ~isempty(fxnname)
    defaultAnswer{2} = fxnname;
    defaultAnswer{3} = sprintf('%s ...', upper(fxnname));
end

% get input
description = inputdlg(inPrompt, inTitle, numLines, defaultAnswer);

% check for user abort
if isempty(description)
    error('description cancelled by user');
end

% read description. Functionname, description and synopsis are required
fxnname = description{2};
if isempty(fxnname) || isempty(description{2})...
        || isempty(description{3}) || isempty(description{4}) ||...
        any(findstr(description{3},'...')) || isempty(regexpi(description{4},fxnname)) 
    error('Username, function name, description and synopsis are required inputs!')
end


% check whether function already exists
fxnCheck = which(fxnname);
if ~isempty(fxnCheck)
    error('%s already exists: %s',fxnname,fxnCheck)
end

% read other input
username = description{1};
desc = description{3};
synopsis = description{4};
% if there are line breaks in inputtext and outputtext: make sure that
% these lines will still be commented!
inputtext = description{5};
if strfind(inputtext,'\n')
    inputtext = regexprep(inputtext,'\\n','\\n%%\\t\\t');
end
outputtext = description{6};
if strfind(outputtext,'\n')
    outputtext = regexprep(outputtext,'\\n','\\n%%\\t\\t\\t');
end

% ask for directory to save the function
dirName = uigetdir(pwd,'select save directory');
if dirName == 0
    error('directory selection cancelled by user')
end
% end of data-gathering stage

%=================

%=================
% WRITE FUNCTION
%=================

% create the filename based on the function name
fsuffix = '.m';
filename = fullfile(dirName,[fxnname fsuffix]);
% end filename creation

% this switch allows users to set their own design preferences.
% Just add case 'myusername'
switch username
    otherwise
        % beginning of file-printing stage
        fid = fopen(filename,'wt');
        fprintf(fid,'function %s\n',synopsis);
        fprintf(fid,'%%%s\n',desc);
        fprintf(fid,'%%\n');
        fprintf(fid,'%% SYNOPSIS: %s\n',synopsis);
        fprintf(fid,'%%\n');
        fprintf(fid,['%% INPUT ',inputtext,'\n']);
        fprintf(fid,'%%\n');
        fprintf(fid,['%% OUTPUT ',outputtext,'\n']);
        fprintf(fid,'%%\n');
        fprintf(fid,'%% REMARKS\n');
        fprintf(fid,'%%\n');
        fprintf(fid,'%% created with MATLAB ver.: %s on %s\n',vers,os);
        fprintf(fid,'%%\n');
        fprintf(fid,'%% created by: %s\n',username);
        fprintf(fid,'%% DATE: %s\n',datetoday);
        fprintf(fid,'%%\n');
        fprintf(fid,'%s\n',repmat('%',1,75));
        fclose(fid);
        % end of file-printing stage
end

% pop up the newly generated file
edit(filename);