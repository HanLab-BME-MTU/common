function varargout = signalProcessingProcessGUI(varargin)
% signalProcessingProcessGUI M-file for signalProcessingProcessGUI.fig
%      signalProcessingProcessGUI, by itself, creates a new signalProcessingProcessGUI or raises the existing
%      singleton*.
%
%      H = signalProcessingProcessGUI returns the handle to a new signalProcessingProcessGUI or the handle to
%      the existing singleton*.
%
%      signalProcessingProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in signalProcessingProcessGUI.M with the given input arguments.
%
%      signalProcessingProcessGUI('Property','Value',...) creates a new signalProcessingProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before signalProcessingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to signalProcessingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help signalProcessingProcessGUI

% Last Modified by GUIDE v2.5 02-Mar-2012 22:38:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @signalProcessingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @signalProcessingProcessGUI_OutputFcn, ...
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


% --- Executes just before signalProcessingProcessGUI is made visible.
function signalProcessingProcessGUI_OpeningFcn(hObject,eventdata,handles,varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});

userData=get(handles.figure1,'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableMovies,'String',userData.ML.movieDataFile_, ...
    'UserData',1:numel(userData.ML.movies_));

movieIndex = funParams.MovieIndex;

% Find any parent process
parentProc = userData.crtPackage.getParent(userData.procID);
if isempty(userData.crtPackage.processes_{userData.procID}) && ~isempty(parentProc)
    % Check existence of all parent processes
    emptyParentProc = any(cellfun(@isempty,userData.crtPackage.processes_(parentProc)));
    if ~emptyParentProc
        % Intersect movie index with channel index of parent processes
        for i = parentProc
            parentParams = userData.crtPackage.processes_{i}.funParams_;
            movieIndex = intersect(movieIndex,parentParams.MovieIndex);
        end
    end
end

if ~isempty(movieIndex)
    movieString = userData.ML.movieDataFile_(movieIndex);
else
    movieString = {};
end

set(handles.listbox_selectedMovies,'String',movieString,'UserData',movieIndex);

% Set tools menu
userData.tools = SignalProcessingProcess.getTools;
nTools=numel(userData.tools);
createToolCheckbox= @(i,type,value) uicontrol(handles.uipanel_tools,...
    'Style','checkbox','Tag',['checkbox_tool' num2str(i)],'Value',0,...
    'Position',[40 10+20*(nTools-i) 250 20],'String',[' ' userData.tools(i).name],'HorizontalAlignment','left');
createToolSettingsButton= @(i,settingsGUI) uicontrol(handles.uipanel_tools,...
    'Style','pushbutton','String','Settings',...
    'Position',[350 10+20*(nTools-i) 100 20],'Tag',['settings_tool' num2str(i)],...
    'Callback',@(h,event) settingsGUI(h,event,guidata(h),i));

resizeGUI(handles,20*nTools);

for i =nTools:-1:1
    createToolCheckbox(i,userData.tools(i).type);
    createToolSettingsButton(i,userData.tools(userData.tools(i).type).GUI);
end

tools= userData.crtProc.funParams_.tools;
for i =1: numel(tools)
    h = findobj(handles.figure1,'Tag',['checkbox_tool' num2str(tools(i).type)]);

    set(h,'Value',1);
    userData.tools(tools(i).type).parameters = tools(i).parameters;    
end
% guidata(handles.figure1,guihandles(handles.figure1));

% Choose default command line output for signalProcessingProcessGUI
handles.output = hObject;

% Update user data and GUI data
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = signalProcessingProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if ishandle(userData.helpFig), delete(userData.helpFig); end
if ishandle(userData.previewFig), delete(userData.previewFig); end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on button press in checkbox_allMovies.
function checkbox_all_Callback(hObject, eventdata, handles)

% Identify listbox and retrieve handles
tokens = regexp(get(hObject,'Tag'),'^checkbox_all(.*)$','tokens');
listbox_available= handles.(['listbox_available' tokens{1}{1}]);
listbox_selected= handles.(['listbox_selected' tokens{1}{1}]);

% Retrieve available properties
availableProps = get(listbox_available, {'String','UserData'});
if isempty(availableProps{1}), return; end

if get(hObject,'Value')
    set(listbox_selected, 'String', availableProps{1},'UserData',availableProps{2});
else
    set(listbox_selected, 'String', {}, 'UserData',[], 'Value',1);
end

% --- Executes on button press in pushbutton_selectMovies.
function pushbutton_select_Callback(hObject, eventdata, handles)

% Identify listbox and retrieve handles
tokens = regexp(get(hObject,'Tag'),'^pushbutton_select(.*)$','tokens');
listbox_available= handles.(['listbox_available' tokens{1}{1}]);
listbox_selected= handles.(['listbox_selected' tokens{1}{1}]);

% Get handles properties
availableProps = get(listbox_available, {'String','UserData'});
selectedProps = get(listbox_selected, {'String','UserData'});
ID = get(listbox_available, 'Value');

% Update selected listbox properties
newChanID = ID(~ismember(availableProps{1}(ID),selectedProps{1}));
selectedString = vertcat(selectedProps{1},availableProps{1}(newChanID));
selectedData = horzcat(selectedProps{2}, availableProps{2}(newChanID));

set(listbox_selected, 'String', selectedString, 'Userdata', selectedData);


% --- Executes on button press in pushbutton_deleteMovies.
function pushbutton_delete_Callback(hObject, eventdata, handles)

% Identify listbox and retrieve handles
tokens = regexp(get(hObject,'Tag'),'^pushbutton_delete(.*)$','tokens');
listbox_selected= handles.(['listbox_selected' tokens{1}{1}]);

% Get selected properties and returin if empty
selectedProps = get(listbox_selected, {'String','UserData','Value'});
if isempty(selectedProps{1}) || isempty(selectedProps{3}),return; end

% Delete selected item
selectedProps{1}(selectedProps{3}) = [ ];
selectedProps{2}(selectedProps{3}) = [ ];
set(listbox_selected, 'String', selectedProps{1},'UserData',selectedProps{2},...
    'Value',max(1,min(selectedProps{3},numel(selectedProps{1}))));


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

% Check user input
if isempty(get(handles.listbox_selectedMovies, 'String'))
    errordlg('Please select at least one input process from ''Available Movies''.','Setting Error','modal')
    return;
end
funParams.MovieIndex = get(handles.listbox_selectedMovies, 'UserData');
    
h = findobj(handles.figure1,'-regexp','Tag','checkbox_tool*');
selectedTools = logical(arrayfun(@(x) get(x,'Value'),h));
userData=get(handles.figure1,'UserData');
funParams.tools=userData.tools(selectedTools);

% Set parameters
processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


function resizeGUI(handles,dh)

toolsUicontrols = get(handles.uipanel_tools,'Children');
resizableHandles=[handles.uipanel_tools,handles.figure1];
arrayfun(@(x) set(x,'Position',get(x,'Position')+[0 0 0 dh]),resizableHandles);
arrayfun(@(x) set(x,'Position',get(x,'Position')+[0 dh 0 0]),toolsUicontrols);

% Fix GUI position/sizebleInput,'Position',get(handles.uipanel_samplableInput,'Position')+[0 0 0 dh]);
uicontrols = {'uipanel_input','axes_help','text_processName','text_copyright'};
for i=1:numel(uicontrols)
    set(handles.(uicontrols{i}),'Position',get(handles.(uicontrols{i}),'Position')+[0 dh 0 0]);
end
