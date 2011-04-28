function varargout = packageGUI(varargin)
% PACKAGEGUI M-file for packageGUI.fig
%      PACKAGEGUI, by itself, creates a new PACKAGEGUI or raises the existing
%      singleton*.
%
%      H = PACKAGEGUI returns the handle to a new PACKAGEGUI or the handle to
%      the existing singleton*.
%
%      PACKAGEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PACKAGEGUI.M with the given input arguments.
%
%      PACKAGEGUI('Property','Value',...) creates a new PACKAGEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before packageGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to packageGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help packageGUI

% Last Modified by GUIDE v2.5 29-Mar-2011 11:31:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @packageGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @packageGUI_OutputFcn, ...
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


% --- Executes just before packageGUI is made visible.
function packageGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% packageGUI(MD)   MD: MovieData object
%
% Useful tools
%
% User Data:
%
%       userData.MD - array of MovieData object
%       userData.package - array of package (same length with userData.MD)
%       userData.crtPackage - the package of current MD
%       userData.id - the id of current MD on board
%
%       userData.dependM - dependency matrix
%       userdata.statusM - GUI status matrix
%       userData.optProcID - optional process ID
%
%       userData.applytoall - array of boolean
%
%       userData.passIconData - pass icon image data
%       userData.errorIconData - error icon image data
%       userData.warnIconData - warning icon image data
%       userData.questIconData - help icon image data
%       userData.colormap - color map
%
%       userData.setFig - array of handles of (multiple) setting figures (may not exist)
%       userData.resultFig - array of handles of (multiple) result figures (may not exist)
%       userData.packageHelpFig - handle of (single) help figure (may not exist)
%       userData.iconHelpFig - handle of (single) help figures (may not exist)
%       userData.processHelpFig - handle of (multiple) help figures (may not exist) 
%       
%
% NOTE:
%   
%   userData.statusM - 1 x m stucture array, m is the number of Movie Data 
%                      this user data is used to save the status of movies
%                      when GUI is switching between different movie(s)
%                   
%   	fields: IconType - the type of status icons, 'pass', 'warn', 'error'
%               Msg - the message displayed when clicking status icons
%               Checked - 1 x n logical array, n is the number of processes
%                         used to save value of check box of each process
%               Visited - logical true or false, if the movie has been
%                         loaded to GUI before 
%

% Load movie data and recycle processes
userfcn_iniPackageGUI;

% --- Outputs from this function are returned to the command line.
function varargout = packageGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% In case the package GUI has been called without argument
userData = get(handles.figure1, 'UserData');
if (isfield(userData,'startMovieSelectorGUI') && userData.startMovieSelectorGUI)
    menu_file_open_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(~, ~, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% userData = get(handles.figure1, 'UserData');
% for i = 1: length(userData.MD)
%     userData.MD(i).save;
% end
delete(handles.figure1);

% --- Executes on button press in pushbutton_status.
function pushbutton_status_Callback(~, ~, handles)
% hObject    handle to pushbutton_status (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

% if movieDataGUI exist
if isfield(userData, 'overviewFig') && ishandle(userData.overviewFig)
    delete(userData.overviewFig)
end

userData.overviewFig = movieDataGUI(userData.MD(userData.id));
set(handles.figure1, 'UserData', userData);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(~, ~, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

for i = 1: length(userData.MD)
    userData.MD(i).save;
end

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')


function switchMovie_Callback(hObject, ~, handles)

userData = get(handles.figure1, 'UserData');
nMovies = length(userData.MD);

switch get(hObject,'Tag')
    case 'pushbutton_left'
        newMovieId = userData.id - 1;
    case 'pushbutton_right'
        newMovieId = userData.id + 1;
    case 'popupmenu_movie'
        newMovieId = get(hObject, 'Value');
    otherwise
end

if (newMovieId==userData.id), return; end

% Save previous movie checkboxes
userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

% Set up new movie GUI parameters
userData.id = mod(newMovieId-1,nMovies)+1;
userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)
set(handles.popupmenu_movie, 'Value', userData.id)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% userfcn_pushbutton_run_common

userfcn_pushbutton_run_common

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
userData = get(handles.figure1,'Userdata');
if isfield(userData, 'MD')
    MD = userData.MD;
else
    delete(handles.figure1);
    return;
end

user_response = questdlg('Do you want to save the current progress?', ...
    'Package Control Panel');
switch lower(user_response)
    case 'yes'
        for i = 1: length(userData.MD)
            userData.MD(i).save;
        end
        delete(handles.figure1);
    case 'no'
        delete(handles.figure1);
    case 'cancel'
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

% Find all figures stored in userData and delete them
if isempty(userData), return; end
userDataFields=fieldnames(userData);
isFig = ~cellfun(@isempty,regexp(userDataFields,'Fig$'));
userDataFigs = userDataFields(isFig);
for i=1:numel(userDataFigs)
     figHandles = userData.(userDataFigs{i});
     validFigHandles = figHandles(ishandle(figHandles)&logical(figHandles));     
     delete(validFigHandles);
end

% msgboxGUI used for error reports
if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
   delete(userData.msgboxGUI) 
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
if strcmp(eventdata.Key, 'leftarrow')
    switchMovie_Callback(handles.pushbutton_left, [], handles);
end
if strcmp(eventdata.Key, 'rightarrow')
    switchMovie_Callback(handles.pushbutton_right, [], handles);
end


% --------------------------------------------------------------------
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

function menu_file_open_Callback(~, ~, handles)
% Call back function of 'New' in menu bar
userData = get(handles.figure1,'Userdata');
if isfield(userData,'MD')
    for i = 1: length(userData.MD)
        userData.MD(i).save;
    end
end
movieSelectorGUI(userData.packageName);
delete(handles.figure1)


% --------------------------------------------------------------------
function menu_file_save_Callback(~, ~, handles)

userData = get(handles.figure1, 'UserData');
userData.MD(userData.id).save;

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')


% --------------------------------------------------------------------
function menu_file_exit_Callback(~, ~, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --------------------------------------------------------------------
function menu_tools_Callback(hObject, eventdata)

handles =guidata(hObject);
userData = get(handles.figure1, 'UserData');
prop=get(hObject,'Tag');
toolID = str2double(prop(length('menu_tools_')+1:end));

toolHandle=userData.crtPackage.tools_(toolID).funHandle;
userData.toolFig(toolID) = toolHandle('mainFig',handles.figure1);

set(handles.figure1, 'UserData', userData);

% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
props=get(hObject,{'Value','Tag'});
procStatus=props{1};
procID = str2double(props{2}(length('pushbutton_clear_')+1:end));

for x = 1: length(userData.MD)
    userData.MD(x).deleteProcess(userData.package(x).processes_{procID})
end


set(handles.figure1, 'UserData', userData);

userfcn_checkAllMovies(procID, procStatus, handles);
userfcn_lampSwitch(procID, procStatus, handles);

% --- Executes on button press in pushbutton_show.
function pushbutton_show_Callback(hObject, ~, handles)

userData = get(handles.figure1, 'UserData');
prop=get(hObject,'Tag');
procID = str2double(prop(length('pushbutton_show_')+1:end));

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    delete(userData.resultFig)
end


% Super-lame way to call the different resultDisplayGUI inputs
% Should work for the moment!
% Modifications should be added to the resultDisplay methods (should be
% generic!!!!)
if isa(userData.crtPackage,'UTrackPackage')
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay(handles.figure1,procID);
else
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay();
end
    
set(handles.figure1, 'UserData', userData);


% --- Executes on button press in pushbutton_set.
function pushbutton_set_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_set (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
prop=get(hObject,'Tag');
procID = str2double(prop(length('pushbutton_set_')+1:end));

%Guess associated proces GUI from process name
crtProc=userData.crtPackage.processClassNames_{procID};
crtProcGUI=str2func([regexprep(crtProc,'(\<[A-Z])','${lower($1)}') 'GUI']);

userData.setFig(procID) = crtProcGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on button press in checkbox.
function checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox
props=get(hObject,{'Value','Tag'});
procStatus=props{1};
procID = str2double(props{2}(length('checkbox_')+1:end));

userfcn_checkAllMovies(procID, procStatus, handles);
userfcn_lampSwitch(procID, procStatus, handles);
