function varargout = detectionProcessGUI(varargin)
%DETECTIONPROCESSGUI M-file for detectionProcessGUI.fig
%      DETECTIONPROCESSGUI, by itself, creates a new DETECTIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = DETECTIONPROCESSGUI returns the handle to a new DETECTIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      DETECTIONPROCESSGUI('Property','Value',...) creates a new DETECTIONPROCESSGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to detectionProcessGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DETECTIONPROCESSGUI('CALLBACK') and DETECTIONPROCESSGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DETECTIONPROCESSGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detectionProcessGUI

% Last Modified by GUIDE v2.5 13-Dec-2011 13:55:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @detectionProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before detectionProcessGUI is made visible.
function detectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:})

% Get current package and process
userData = get(handles.figure1, 'UserData');

% Get current process constructer, set-up GUIs and mask refinement process
% constructor
     
userData.subProcClassNames = eval([userData.crtProcClassName '.getConcreteClasses']);
validClasses = cellfun(@(x)exist(x,'class')==8,userData.subProcClassNames);
userData.subProcClassNames = userData.subProcClassNames(validClasses);
userData.subProcConstr = cellfun(@(x) str2func(x),userData.subProcClassNames,'Unif',0);
userData.subProcGUI = cellfun(@(x) eval([x '.GUI']),userData.subProcClassNames,'Unif',0);
subProcNames = cellfun(@(x) eval([x '.getName']),userData.subProcClassNames,'Unif',0);
popupMenuProcName = vertcat(subProcNames,{'Choose a detection method'});

% Set up input channel list box
if isempty(userData.crtProc)
    value = numel(userData.subProcClassNames)+1;
    set(handles.pushbutton_set, 'Enable', 'off');
else
    value = find(strcmp(userData.crtProc.getName,subProcNames));
end

existSubProc = @(proc) any(cellfun(@(x) isa(x,proc),userData.MD.processes_));
for i=find(cellfun(existSubProc,userData.subProcClassNames'))
  popupMenuProcName{i} = ['<html><b>' popupMenuProcName{i} '</b></html>'];
end

set(handles.popupmenu_methods, 'String', popupMenuProcName,...
    'Value',value)

% Choose default command line output for detectionProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = detectionProcessGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
delete(handles.figure1);


% --- Executes on selection change in popupmenu_methods.
function popupmenu_methods_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_methods contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_methods
content = get(hObject, 'string');
if get(hObject, 'Value') == length(content)
    set(handles.pushbutton_set, 'Enable', 'off')
else
    set(handles.pushbutton_set, 'Enable', 'on')
end

% --- Executes on button press in pushbutton_set.
function pushbutton_set_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
segProcID = get(handles.popupmenu_methods, 'Value');
subProcGUI =userData.subProcGUI{segProcID};
subProcGUI('mainFig',userData.mainFig,userData.procID,...
    'procConstr',userData.subProcConstr{segProcID},...
    'procClassName',userData.subProcClassNames{segProcID});
delete(handles.figure1);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
