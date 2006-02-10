function gencodetemplate
% GENCODETEMPLATE creates a template for new functions
%
% DESCRIPTION: This function generates a template based on user-defined
% variables for new functions.  The point is to improve the documentation
% process during the software development phase.
%
%
%
% SYNOPSIS: gencodetemplate
%
%
%
% INPUT: Input is derived via user-defined variables .
%
%
%
% OUTPUT: The output is an m-file function filled in with user-defined
% variables.
%
%
%
% MATLAB VERSION (originally written on): 7.1.0.246 (R14) Service Pack 3
%
%
%
% USERNAME: jkunken
%
%
%

disp(sprintf('LCCB code template generator - modified by jonas\n'));

% data-gathering stage begins below
% getenv() and version are used to derive environmental variables

% find username. On Linux, there is no environment variable, therefore
% we'll have to ask, anyway.
username = getenv('username');
if isempty(username)
    % this should only happen on Linux, but better be safe.
    usernameFeedback = 'y';
else
    disp(sprintf('Your username appears to be %s.\n',username));
    usernameFeedback = input('Is this correct? (y/n)\n','s');
end
% allow user to change name
if strcmp(usernameFeedback,'y')
    % all is well
else
    username = input('Please enter your username:\n','s');
end
% check that we do have a non-empty username
if isempty(username)
    error('no username provided')
end

% get version, OS, date
vers    = version;
os      = getenv('OS');
datetoday = date;
% ask for description
inPrompt = {'Function name',...
    'Description: ''FUNCTIONNAME does ...'' (captitalized function name)',...
    'Synopsis: ''[output1, output2] = functionname(input1, input2)',...
    'Description of input arguments (use \n for line breaks)',...
    'Description of output arguments (use \n for line breaks)'};
inTitle = 'Please describe your function';
numLines = repmat([1,100],5,1);
[defaultAnswer{1:5}] = deal('');
description = inputdlg(inPrompt, inTitle, numLines, defaultAnswer);
% check for user abort
if isempty(description)
    error('description cancelled by user');
end
% read description. Functionname, description and synopsis are required
if isempty(description{1}) || isempty(description{2}) || isempty(description{3})
    error('Function name, description and synopsis are required inputs!')
end

fxnname = description{1};
desc = description{2};
synopsis = description{3};
% if there are line breaks in inputtext and outputtext: make sure that
% these lines will still be commented!
inputtext = description{4};
if strfind(inputtext,'\n')
    inputtext = regexprep(inputtext,'\\n','\\n%%\\t\\t');
end
outputtext = description{5};
if strfind(outputtext,'\n')
    outputtext = regexprep(outputtext,'\\n','\\n%%\\t\\t\\t');
end
dirName = uigetdir(pwd,'select save directory');
if dirName == 0
    error('directory selection cancelled by user')
end
% end of data-gathering stage

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
        fprintf(fid,'%%\n');
        fclose(fid);
        % end of file-printing stage
end

disp(sprintf('Your code template is now ready\n'));

% pop up the newly generated file
edit(filename);

%clean up all of the mess we made in the workspace
clear;