function varargout = segmentationProcessGUI(varargin)
%SEGMENTATIONPROCESSGUI M-file for segmentationProcessGUI.fig
%      SEGMENTATIONPROCESSGUI, by itself, creates a new SEGMENTATIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = SEGMENTATIONPROCESSGUI returns the handle to a new SEGMENTATIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      SEGMENTATIONPROCESSGUI('Property','Value',...) creates a new SEGMENTATIONPROCESSGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to segmentationProcessGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SEGMENTATIONPROCESSGUI('CALLBACK') and SEGMENTATIONPROCESSGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SEGMENTATIONPROCESSGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segmentationProcessGUI

% Last Modified by GUIDE v2.5 28-Sep-2011 16:04:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segmentationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @segmentationProcessGUI_OutputFcn, ...
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


% --- Executes just before segmentationProcessGUI is made visible.
function segmentationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.setFig(procID) = segmentationProcessGUI('mainFig',handles.figure1, procID);
%
% Available tools 
% UserData data:
%       userData.MD - 1x1 the current movie data
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main figure
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%
%
%       userData.segProc - cell array of segmentation processes, created
%                          after setting up segmentation processes
%
%
%       userData.procSetting - cell array of set-up GUIs of available
%                              processes
%       userData.procName - cell array of available segmentation processes
%       userData.procConstr - constructor of current process
%
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:})

% Get current package and process
userData = get(handles.figure1, 'UserData');

% Get current process constructer, set-up GUIs and mask refinement process
% constructor
     
userData.subProcClassNames = eval([func2str(userData.procConstr) '.getConcreteClasses']);
validClasses = cellfun(@(x)exist(x,'class')==8,userData.subProcClassNames);
userData.subProcClassNames = userData.subProcClassNames(validClasses);
userData.subProcConstr = cellfun(@(x) str2func(x),userData.subProcClassNames,'Unif',0);
userData.subProcGUI = cellfun(@(x) eval([x '.GUI']),userData.subProcClassNames,'Unif',0);
subProcNames = cellfun(@(x) eval([x '.getName']),userData.subProcClassNames,'Unif',0);
popupMenuProcName = vertcat(subProcNames,{'Choose a segmentation method'});

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

set(handles.popupmenu_segmentationMethods, 'String', popupMenuProcName,...
    'Value',value)

% Choose default command line output for segmentationProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = segmentationProcessGUI_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on selection change in popupmenu_segmentationMethods.
function popupmenu_segmentationMethods_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_segmentationMethods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_segmentationMethods contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_segmentationMethods
content = get(hObject, 'string');
if get(hObject, 'Value') == length(content)
    set(handles.pushbutton_set, 'Enable', 'off')
else
    set(handles.pushbutton_set, 'Enable', 'on')
end

% --- Executes on button press in pushbutton_set.
function pushbutton_set_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
segProcID = get(handles.popupmenu_segmentationMethods, 'Value');
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
