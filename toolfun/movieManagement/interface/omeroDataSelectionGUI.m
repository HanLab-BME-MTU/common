function varargout = omeroDataSelectionGUI(varargin)
% OMERODATASELECTIONGUI MATLAB code for omeroDataSelectionGUI.fig
%      OMERODATASELECTIONGUI, by itself, creates a new OMERODATASELECTIONGUI or raises the existing
%      singleton*.
%
%      H = OMERODATASELECTIONGUI returns the handle to a new OMERODATASELECTIONGUI or the handle to
%      the existing singleton*.
%
%      OMERODATASELECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OMERODATASELECTIONGUI.M with the given input arguments.
%
%      OMERODATASELECTIONGUI('Property','Value',...) creates a new OMERODATASELECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before omeroDataSelectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to omeroDataSelectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help omeroDataSelectionGUI

% Last Modified by GUIDE v2.5 11-Nov-2013 21:21:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @omeroDataSelectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @omeroDataSelectionGUI_OutputFcn, ...
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


% --- Executes just before omeroDataSelectionGUI is made visible.
function omeroDataSelectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to omeroDataSelectionGUI (see VARARGIN)

set(handles.text_copyright, 'String', getLCCBCopyright())
global session

% List projects
projects = getProjects(session);        
projectNames = arrayfun(@(x) char(x.getName().getValue()), projects,...
    'UniformOutput', false);
projectNames = [{'none'}; projectNames];
set(handles.popupmenu_project, 'Value', 1, 'String', projectNames,...
    'UserData', projects);

% List orphaned datasets
parameters = omero.sys.ParametersI();
parameters.orphan();
parameters.leaves()
datasets = getObjects(session, 'dataset', [], parameters);
datasetNames = arrayfun(@(x) char(x.getName().getValue()), datasets,...
    'UniformOutput', false);
datasetNames = [{'none'}; datasetNames];
set(handles.popupmenu_dataset, 'Value', 1, 'String', datasetNames,...
    'UserData', datasets);
refreshImageList(handles)

% Save orphaned dataset
userData = get(hObject, 'UserData');
userData.orphaned_datasets = datasets;
set(hObject, 'UserData', userData);

% Choose default command line output for omeroDataSelectionGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes omeroDataSelectionGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = omeroDataSelectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function refreshImageList(handles)

% Read image list from selected dataset
props = get(handles.popupmenu_dataset, {'Value', 'UserData'});
if props{1} == 1,
    set(handles.listbox_images, 'Value', 1, 'String', '');
    images = [];
else
    images = toMatlabList(props{2}(props{1} - 1).linkedImageList);
end

% Update image list
imageNames = arrayfun(@(x) char(x.getName().getValue()), images,...
    'UniformOutput', false);
set(handles.listbox_images, 'Value', 1, 'String', imageNames,...
    'UserData', images);

% --- Executes on selection change in popupmenu_project.
function popupmenu_project_Callback(hObject, eventdata, handles)

% Read dataset list from selected project
props = get(hObject, {'Value', 'UserData'});
if props{1} == 1,
    userData = get(handles.figure1, 'UserData');
    datasets = userData.orphaned_datasets;
else
    datasets = toMatlabList(props{2}(props{1} - 1).linkedDatasetList);
end

% Update datasets drop down menu list
datasetNames = arrayfun(@(x) char(x.getName().getValue()), datasets,...
    'UniformOutput', false);
datasetNames = [{'none'}; datasetNames];
set(handles.popupmenu_dataset, 'Value', 1, 'String', datasetNames,...
    'UserData', datasets);
refreshImageList(handles)

% --- Executes on selection change in popupmenu_dataset.
function popupmenu_dataset_Callback(hObject, eventdata, handles)
refreshImageList(handles)
