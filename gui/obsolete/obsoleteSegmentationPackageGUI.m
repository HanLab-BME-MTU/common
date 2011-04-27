function varargout = obsoleteSegmentationPackageGUI(varargin)
% OBSOLETESEGMENTATIONPACKAGEGUI M-file for obsoleteSegmentationPackageGUI.fig
%      OBSOLETESEGMENTATIONPACKAGEGUI, by itself, creates a new OBSOLETESEGMENTATIONPACKAGEGUI or raises the existing
%      singleton*.
%
%      H = OBSOLETESEGMENTATIONPACKAGEGUI returns the handle to a new OBSOLETESEGMENTATIONPACKAGEGUI or the handle to
%      the existing singleton*.
%
%      OBSOLETESEGMENTATIONPACKAGEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OBSOLETESEGMENTATIONPACKAGEGUI.M with the given input arguments.
%
%      OBSOLETESEGMENTATIONPACKAGEGUI('Property','Value',...) creates a new OBSOLETESEGMENTATIONPACKAGEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before obsoleteSegmentationPackageGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to obsoleteSegmentationPackageGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help obsoleteSegmentationPackageGUI

% Last Modified by GUIDE v2.5 26-Apr-2011 10:33:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @obsoleteSegmentationPackageGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @obsoleteSegmentationPackageGUI_OutputFcn, ...
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


% --- Executes just before obsoleteSegmentationPackageGUI is made visible.
function obsoleteSegmentationPackageGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% biosensorsPackageGUI(MD)   MD: MovieData object
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
%       userData.passIconData - pass icon image data
%       userData.errorIconData - error icon image data
%       userData.warnIconData - warning icon image data
%       userData.questIconData - help icon image data
%       userData.colormap - color map
%
%       userData.setFig - array of handles of (multiple) setting figures (may not exist)
%       userData.resultFig - array of handles of (multiple) result figures (may not exist)
%       userData.setupMovieDataFig - handle of (single) setupMovieData figure (may not exist)
%       userData.overviewFig - handles of (single) overviewMovieDataGUI figure (may not exist)
%       userData.packageHelpFig - handle of (single) help figure (may not exist)
%       userData.iconHelpFig - handle of (single) help figures (may not exist)
%       userData.processHelpFig - handle of (multiple) help figures (may not exist) 
%       userData.msgboxGUI - handle of message box
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
userfcn_iniPackage_commonCode;



% --- Outputs from this function are returned to the command line.
function varargout = obsoleteSegmentationPackageGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% In case the GUI has been called without argument
userData = get(handles.figure1,'Userdata');
if (isfield(userData,'startMovieSelectorGUI') && userData.startMovieSelectorGUI)
    menu_file_open_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
for i = 1: length(userData.MD)
    userData.MD(i).save
end
delete(handles.figure1);


% --- Executes on button press in pushbutton_status.
function pushbutton_status_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

% if movieDataGUI exist
if isfield(userData, 'overviewFig') && ishandle(userData.overviewFig)
    delete(userData.overviewFig)
end

userData.overviewFig = movieDataGUI('mainFig',handles.figure1, 'overview', userData.MD(userData.id));
set(handles.figure1, 'UserData', userData);

% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

for i = 1: length(userData.MD)
    userData.MD(i).save
end

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')


% --- Executes on button press in checkbox_1.
function checkbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_1
userfcn_checkAllMovies(1, get(hObject,'value'), handles);
userfcn_lampSwitch(1, get(hObject,'value'), handles);

% --- Executes on button press in pushbutton_set_1.
function pushbutton_set_1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_set_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
procID = 1;
userData.setFig(procID) = segmentationProcessGUI('mainFig',handles.figure1,procID);
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_show_1.
function pushbutton_show_1_Callback(hObject, eventdata, handles)

procID = 1;
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'resultFig') && ishandle(userData.resultFig)
    
    delete(userData.resultFig)
end
%     userData.resultFig = userData.crtPackage.processes_{procID}.showResult;
    userData.resultFig = userData.crtPackage.processes_{procID}.resultDisplay;


set(handles.figure1, 'UserData', userData);


% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)

userfcn_pushbutton_run_common



% --- Executes on button press in checkbox_forcerun.
function checkbox_forcerun_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_forcerun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_forcerun


% --- Executes on button press in pushbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)
% userData.id
% userData.crtPackage
%
userData = get(handles.figure1, 'UserData');
l = length(userData.MD);

userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = userData.id - 1;

if userData.id < 1
   userData.id = l;
end

userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up movie explorer
set(handles.popupmenu_movie, 'Value', userData.id)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end


% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
l = length(userData.MD);

userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = userData.id + 1;

if userData.id > l
   userData.id = mod(userData.id, l);
end

userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up movie explorer
set(handles.popupmenu_movie, 'Value', userData.id)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end


% --- Executes on selection change in popupmenu_movie.
function popupmenu_movie_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if get(hObject, 'Value') == userData.id
   return 
end

l = length(userData.MD);
userData.statusM(userData.id).Checked = userfcn_saveCheckbox(handles);

userData.id = get(hObject, 'Value');
userData.crtPackage = userData.package(userData.id);
set(handles.figure1, 'UserData', userData)

% Set up GUI
if userData.statusM(userData.id).Visited
   userfcn_updateGUI(handles, 'refresh') 
else
   userfcn_updateGUI(handles, 'initialize') 
end


% --- Executes during object creation, after setting all properties.
function popupmenu_movie_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_runall.
function checkbox_runall_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_runall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_runall


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
if strcmp(eventdata.Key, 'leftarrow')
    pushbutton_left_Callback(handles.pushbutton_left, [], handles);
end
if strcmp(eventdata.Key, 'rightarrow')
    pushbutton_right_Callback(handles.pushbutton_right, [], handles);
end


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
    'Segmentation Package Control Panel');
switch lower(user_response)
    case 'yes'
        for i = 1: length(userData.MD)
            userData.MD(i).save
        end
        delete(handles.figure1);
    case 'no'
        delete(handles.figure1);
    case 'cancel'
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
% setFlag = getappdata(hObject, 'setFlag');

% Delete setting figures (multiple)
if isfield(userData, 'setFig')
    for i = 1: length(userData.setFig)
        if userData.setFig(i)~=0 && ishandle(userData.setFig(i))
            delete(userData.setFig(i))
        end
    end
end

% Close result figures (multiple)
if isfield(userData, 'resultFig') && userData.resultFig~=0 && ishandle(userData.resultFig)

	delete(userData.resultFig(i))
end

% If open, delete MovieData Overview GUI figure (single)
if isfield(userData, 'setupMovieDataFig') && ishandle(userData.setupMovieDataFig)
   delete(userData.setupMovieDataFig) 
end

% If open, delete MovieData Overview GUI figure (single)
if isfield(userData, 'overviewFig') && ishandle(userData.overviewFig)
   delete(userData.overviewFig) 
end

% Delete pre-defined package help dialog (single)
if isfield(userData, 'packageHelpFig') && ishandle(userData.packageHelpFig)
   delete(userData.packageHelpFig) 
end

% Delete pre-defined icon help dialog (single)
if isfield(userData, 'iconHelpFig') && ishandle(userData.iconHelpFig)
   delete(userData.iconHelpFig) 
end

% Delete pre-defined process help dialogssetting figures (multiple)
if isfield(userData, 'processHelpFig')
    for i = 1: length(userData.processHelpFig)
        if userData.processHelpFig(i)~=0 && ishandle(userData.processHelpFig(i))
            delete(userData.processHelpFig(i))
        end
    end
end

% msgboxGUI used for error reports
if isfield(userData, 'msgboxGUI') && ishandle(userData.msgboxGUI)
   delete(userData.msgboxGUI) 
end


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_lccb_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_lccb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_about_update_Callback(hObject, eventdata, handles)
status = web('http://lccb.hms.harvard.edu/software.html', '-browser');
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
function menu_about_lccb_Callback(hObject, eventdata, handles)
% hObject    handle to menu_about_lccb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
status = web('http://lccb.hms.harvard.edu/', '-browser');
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
function menu_file_open_Callback(hObject, eventdata, handles)
% Call back function of 'New' in menu bar
userData = get(handles.figure1,'Userdata');
if isfield(userData,'MD')
    for i = 1: length(userData.MD)
        userData.MD(i).save
    end
end
movieSelectorGUI(userData.packageName);
delete(handles.figure1)


% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');
userData.MD(userData.id).save

set(handles.text_body3, 'Visible', 'on')
pause(1)
set(handles.text_body3, 'Visible', 'off')


% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);


% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)
% userData = get(handles.figure1, 'UserData');
% l = size(userData.dependM, 1);

% switch get(hObject,'value')
%     case 1
%         userfcn_enable(1:length(userData.dependM(:,1)),'on',handles,true);
        
        % track checked status
%         for x = 1: length(userData.MD)
%             userData.statusM(x).Checked = ones(1, l);
%         end
        
%     case 0
%         k = [];
%         for i = 1: size(userData.dependM, 1)
%             if ~isempty(userData.crtPackage.processes_{i}) && ...
%                 userData.crtPackage.processes_{i}.success_ 
%                 k = [k, i];
%             end
%             eval( ['set(handles.checkbox_',num2str(i),',''value'',0)']  );
%         end
%         tempDependM = userData.dependM;
%         tempDependM(:,k) = zeros(size(userData.dependM,1),length(k));
%         userfcn_enable(find(any(tempDependM,2)),'off',handles,true);
        
        % track checked status
%         for x = 1: length(userData.MD)
%             userData.statusM(x).Checked = zeros(1, l);
%         end        
        
% end

% set(handles.figure1, 'UserData', userData)



function edit_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path as text
%        str2double(get(hObject,'String')) returns contents of edit_path as a double


% --- Executes during object creation, after setting all properties.
function edit_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
