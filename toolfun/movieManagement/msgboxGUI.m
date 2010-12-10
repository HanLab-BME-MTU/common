function varargout = msgboxGUI(varargin)
% MSGBOXGUI M-file for msgboxGUI.fig
%      MSGBOXGUI, by itself, creates a new MSGBOXGUI or raises the existing
%      singleton*.
%
%      H = MSGBOXGUI returns the handle to a new MSGBOXGUI or the handle to
%      the existing singleton*.
%
%      MSGBOXGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSGBOXGUI.M with the given input arguments.
%
%      MSGBOXGUI('Property','Value',...) creates a new MSGBOXGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before msgboxGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to msgboxGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help msgboxGUI

% Last Modified by GUIDE v2.5 08-Jun-2010 13:35:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @msgboxGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @msgboxGUI_OutputFcn, ...
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


% --- Executes just before msgboxGUI is made visible.
function msgboxGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% fhandle = msgboxGUI
%
% fhandle = msgboxGUI(paramName, 'paramValue', ... )
%
% Help dialog GUI
% 
% Parameter Field Names:
%       'text' -> Help text
%       'title' -> Title of text box
%       'name'-> Name of dialog box
%

userData = get(handles.figure1, 'UserData');
% Choose default command line output for msgboxGUI
handles.output = hObject;

if nargin > 3
    for i = 1:2:(nargin-3)
        switch lower(varargin{i})
            case 'text'
                set(handles.edit_1, 'String', varargin{i+1})
            case 'name'
                set(hObject, 'Name', varargin{i+1})
            case 'title'
                set(handles.text_title, 'string', varargin{i+1})
            otherwise
                error('User-defined: error calling msgboxGUI. Refer to msgboxGUI.m')
        end
    end
end

% Update handles structure
set(handles.figure1, 'Userdata',userData)
uicontrol(handles.pushbutton_done)
guidata(hObject, handles);

% UIWAIT makes msgboxGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = msgboxGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_1 as text
%        str2double(get(hObject,'String')) returns contents of edit_1 as a double


% --- Executes during object creation, after setting all properties.
function edit_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
delete(handles.figure1)


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    delete(hObject)
end
