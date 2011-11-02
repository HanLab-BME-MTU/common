function varargout = movieSelectorGUI(varargin)
% MOVIESELECTORGUI M-file for movieSelectorGUI.fig
%      MOVIESELECTORGUI, by itself, creates a new MOVIESELECTORGUI or raises the existing
%      singleton*.
%
%      H = MOVIESELECTORGUI returns the handle to a new MOVIESELECTORGUI or the handle to
%      the existing singleton*.
%
%      MOVIESELECTORGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVIESELECTORGUI.M with the given input arguments.
%
%      MOVIESELECTORGUI('Property','Value',...) creates a new MOVIESELECTORGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before movieSelectorGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to movieSelectorGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help movieSelectorGUI

% Last Modified by GUIDE v2.5 22-Sep-2011 15:27:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @movieSelectorGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @movieSelectorGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before movieSelectorGUI is made visible.
function movieSelectorGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% Useful tools:
% 
% User Data:
%
%   userData.MD - new or loaded MovieData object
%   userData.ML - newly saved or loaded MovieList object
%
%   userData.userDir - default open directory
%   userData.colormap - color map (used for )
%   userData.questIconData - image data of question icon
%
%   userData.packageGUI - the name of package GUI
%
%   userData.newFig - handle of new movie set-up GUI
%   userData.iconHelpFig - handle of help dialog
%   userData.msgboxGUI - handle of message box GUI

ip = inputParser;
ip.addRequired('hObject',@ishandle);
ip.addRequired('eventdata',@(x) isstruct(x) || isempty(x));
ip.addRequired('handles',@isstruct);
ip.addParamValue('packageName','',@ischar);
ip.addParamValue('MD',[],@(x) isempty(x) || isa(x,'MovieData'));
ip.parse(hObject,eventdata,handles,varargin{:});

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');

% Choose default command line output for setupMovieDataGUI
handles.output = hObject;

% other user data set-up
userData.MD = [ ];
userData.ML = [ ];
userData.userDir = pwd;
userData.newFig=-1;
userData.msgboxGUI=-1;
userData.iconHelpFig =-1;

% Check packages availability
packageRadioButtons  = get(handles.uipanel_package,'Children');
packageList = get(packageRadioButtons,'UserData');
isValidPackage=logical(cellfun(@(x) exist(x,'class'),packageList));
if isempty(isValidPackage), 
    warndlg('No package found! Please make sure you properly added the installation directory to the path (see user''s manual).',...
        'Movie Selector','modal'); 
end
set(handles.uipanel_package,'SelectionChangeFcn','');

% Grey out radio buttons corresponding to non-available packages
invalidRadioButtons = packageRadioButtons(~isValidPackage);
set(invalidRadioButtons,'Enable','off');
if ~isempty(ip.Results.packageName)
    preSelectedPackage=strcmp(packageList,ip.Results.packageName);
    set(packageRadioButtons(preSelectedPackage),'Value',1.0);
end

% Test a package preselection and update the corresponding radio button
if ~isempty(ip.Results.MD)
    userData.MD=ip.Results.MD;
    nMovies = numel(userData.MD);
    contentlist =  arrayfun(@(x) [x.movieDataPath_ filesep ...
        x.movieDataFileName_],userData.MD,'UniformOutput',false);
    
    % Refresh movie list box in movie selector panel
    set(handles.listbox_movie, 'String', contentlist,'Value',1);
    title = sprintf('Movie List: %g/%g movie(s)',1,nMovies);
    set(handles.text_movie_1, 'String', title)
end

% Load help icon from dialogicons.mat
load lccbGuiIcons.mat
supermap(1,:) = get(hObject,'color');

userData.colormap = supermap;
userData.questIconData = questIconData;

set(handles.figure1,'CurrentAxes',handles.axes_help);
Img = image(questIconData);
set(hObject,'colormap',supermap);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);

if openHelpFile, set(Img, 'UserData', struct('class',mfilename)); end

% Save userdata
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = movieSelectorGUI_OutputFcn(hObject, eventdata, handles) 
% %varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1)

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check a movie and a package are selected
if isempty(get(handles.listbox_movie, 'String'))
   warndlg('Please select at least one movie to continue.', 'Movie Selector', 'modal')
   return
end
if isempty(get(handles.uipanel_package, 'SelectedObject'))
   warndlg('Please select a package to continue.', 'Movie Selector', 'modal')
   return
end

% Retrieve the ID of the selected button and call the appropriate
userData = get(handles.figure1, 'userdata');
selectedPackage=get(get(handles.uipanel_package, 'SelectedObject'),'UserData');
close(handles.figure1);
packageGUI(selectedPackage,userData.MD);


% --- Executes on selection change in listbox_movie.
function listbox_movie_Callback(hObject, eventdata, handles)

contentlist = get(handles.listbox_movie, 'String');
if isempty(contentlist)
    title = sprintf('Movie List: 0/0 movie(s)');
else
    title = sprintf('Movie List: %s/%s movie(s)',...
        num2str(get(handles.listbox_movie, 'Value')),...
        num2str(length(contentlist)));
end
set(handles.text_movie_1, 'String', title)

% --- Executes on button press in pushbutton_new.
function pushbutton_new_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
% if movieDataGUI exist, delete it
if ishandle(userData.newFig), delete(userData.newFig); end
userData.newFig = movieDataGUI('mainFig',handles.figure1);
set(handles.figure1,'UserData',userData);

% --- Executes on button press in pushbutton_prepare.
function pushbutton_prepare_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
% if preparation GUI exist, delete it
if ishandle(userData.newFig), delete(userData.newFig); end
userData.newFig = dataPreparationGUI('mainFig',handles.figure1);
set(handles.figure1,'UserData',userData);

% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'Userdata');
contentlist = get(handles.listbox_movie,'String');
% Return if list is empty
if isempty(contentlist), return;end
num = get(handles.listbox_movie,'Value');

% Delete channel object
delete(userData.MD(num))
userData.MD(num) = [];

% Refresh listbox_channel
contentlist(num) = [ ];
set(handles.listbox_movie,'String',contentlist);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if num>length(contentlist) && num>1
    set(handles.listbox_movie,'Value',length(contentlist));
end
if isempty(contentlist)
    title = sprintf('Movie List: 0/0 movie(s)');
    set(handles.text_movie_1, 'String', title)
else
    title = sprintf('Movie List: %s/%s movie(s)',...
        num2str(get(handles.listbox_movie, 'Value')),...
        num2str(length(contentlist)));
    set(handles.text_movie_1, 'String', title)    
end

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);

% --- Executes on button press in pushbutton_detail.
function pushbutton_detail_Callback(hObject, eventdata, handles)

% Return if no movie 
props=get(handles.listbox_movie, {'String','Value'});
if isempty(props{1}), return; end

userData = get(handles.figure1, 'UserData');
% if movieDataGUI exist, delete it
if ishandle(userData.newFig), delete(userData.newFig); end
userData.newFig = movieDataGUI(userData.MD(props{2}));
set(handles.figure1,'UserData',userData);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
% Delete new window

if ishandle(userData.newFig), delete(userData.newFig); end
if ishandle(userData.iconHelpFig), delete(userData.iconHelpFig); end
if ishandle(userData.msgboxGUI), delete(userData.msgboxGUI); end

% --- Executes on button press in pushbutton_open.
function pushbutton_open_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

[filename, pathname] = uigetfile('*.mat','Select a movie MAT file', ...
    userData.userDir);
if ~any([filename pathname]), return; end
userData.userDir = pathname;

% Check if reselect the movie data that is already in the listbox
contentlist = get(handles.listbox_movie, 'String');

if any(strcmp([pathname filename], contentlist))
    errordlg('This movie has already been selected.','Error','modal');
    return
end

try
    M = MovieObject.load([pathname filename]);
catch ME
    msg = sprintf('Movie: %s\n\nError: %s\n\nMovie is not successfully loaded. Please refer to movie detail and adjust your data.', [pathname filename],ME.message);
    errordlg(msg, 'Movie error','modal');
    return
end

switch class(M)
    case 'MovieData'
        userData.MD = cat(2, userData.MD, M);
        contentlist{end+1} = [pathname filename];
        
    case 'MovieList'        
        % Find duplicate movie data in list box
        movieDataFile = M.movieDataFile_;
        index = 1: length(movieDataFile);
        index = index( ~cellfun(@(z)any(z), cellfun(@(x)strcmp(x, contentlist), movieDataFile, 'UniformOutput', false), 'UniformOutput', true) );
        
        if isempty(index)
            msg = sprintf('All movie data in movie list file %s has already been added to the movie list box.', M.movieListFileName_);
            warndlg(msg,'Warning','modal');
            return
        end
        
        % Reload movie data filenames in case they have been relocated
        % during sanity check        
        movieException = M.sanityCheck();
        movieException = movieException(index);
        errorME = find(~cellfun(@isempty,movieException));
        healthMD = find(cellfun(@isempty,movieException));
        
        % Error movie index
        if ~isempty(errorME)
            filemsg = '';
            for i = errorME
                filemsg = cat(2, filemsg, sprintf('Movie %d:  %s\nError:  %s\n\n', index(i), movieDataFile{index(i)}, movieException{i}.message));
            end
            msg = sprintf('The following movie(s) cannot be sucessfully loaded:\n\n%s', filemsg);
            titlemsg = sprintf('Movie List: %s', [pathname filename]);
            userData.msgboxGUI = msgboxGUI('title',titlemsg,'text', msg);
        end
        
        % Healthy Movie Data
        if ~isempty(healthMD)
            userData.ML = horzcat(userData.ML, M);
            userData.MD = horzcat(userData.MD,M.movies_{index(healthMD)});
            contentlist = horzcat(contentlist',M.movieDataFile_(index(healthMD)));
        end             
    otherwise
        error('User-defined: varable ''type'' does not have an approprate value.')
end

% Refresh movie list box in movie selector panel
set(handles.listbox_movie, 'String', contentlist, 'Value', length(contentlist))
title = sprintf('Movie List: %s/%s movie(s)', num2str(get(handles.listbox_movie, 'Value')), num2str(length(contentlist)));
set(handles.text_movie_1, 'String', title)

set(handles.figure1, 'UserData', userData);

function menu_about_Callback(hObject, eventdata, handles)

status = web(get(hObject,'UserData'), '-browser');
if status
    switch status
        case 1
            msg = 'System default web browser is not found.';
        case 2
            msg = 'System default web browser is found but could not be launched.';
        otherwise
            msg = 'Fail to open browser for unknown reason.';
    end
    warndlg(msg,'Fail to open browser','modal');
end

% --------------------------------------------------------------------
function menu_file_quit_Callback(hObject, eventdata, handles)
delete(handles.figure1)

% --- Executes on button press in pushbutton_deleteall.
function pushbutton_deleteall_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'Userdata');

contentlist = get(handles.listbox_movie,'String');
% Return if list is empty
if isempty(contentlist), return; end
 
user_response = questdlg(['Are you sure to delete all the '...
    num2str(length(contentlist)) ' movie(s) in the listbox?'], ...
    'Movie Listbox', 'Yes','No','Yes');

if strcmpi('no', user_response), return; end

% Delete channel object
userData.MD = [];

% Refresh listbox_channel
set(handles.listbox_movie,'String',{}, 'Value',1);
set(handles.text_movie_1, 'String','Movie List: 0/0 movie(s)')

set(handles.figure1, 'Userdata', userData)
guidata(hObject, handles);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

contentList = get(handles.listbox_movie, 'String');
if isempty(contentList)
    warndlg('No movie selected. Please create new movie data or open existing movie data or movie list.', 'No Movie Selected', 'modal')
    return
end

if isempty(userData.ML)
    movieListPath = [userData.userDir filesep];
    movieListFileName = 'movieDataList.mat';
else
    movieListPath = userData.ML(end).movieListPath_;
    movieListFileName = userData.ML(end).movieListFileName_;
end

% Ask user where to save the movie data file
[filename,path] = uiputfile('*.mat','Find a place to save your movie list',...
             [movieListPath filesep movieListFileName]);
         
if ~any([filename,path]), return; end

% Ask user where to select the output directory of the
outputDir = uigetdir(path,'Select a directory to store the processes output');
if isequal(outputDir,0), return; end

try
    ML = MovieList(contentList, outputDir,'movieListPath_', path,...
        'movieListFileName_', filename);
catch ME
    msg = sprintf('%s\n\nMovie list is not saved.', ME.message);
    errordlg(msg, 'Movie List Error', 'modal')
    return
end

% Save the movie list
ML.save();

% --------------------------------------------------------------------
function menu_tools_crop_Callback(hObject, eventdata, handles)

% Return if no movie 
props=get(handles.listbox_movie, {'String','Value'});
if isempty(props{1}), return; end

userData = get(handles.figure1, 'UserData');
if ishandle(userData.newFig), delete(userData.newFig); end
userData.newFig = cropMovieGUI(userData.MD(props{2}),'mainFig',handles.figure1);
set(handles.figure1,'UserData',userData);