function varargout = channelGUI(varargin)
% CHANNELGUI M-file for channelGUI.fig
%      CHANNELGUI, by itself, creates a new CHANNELGUI or raises the existing
%      singleton*.
%
%      H = CHANNELGUI returns the handle to a new CHANNELGUI or the handle to
%      the existing singleton*.
%
%      CHANNELGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHANNELGUI.M with the given input arguments.
%
%      CHANNELGUI('Property','Value',...) creates a new CHANNELGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before channelGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to channelGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help channelGUI

% Last Modified by GUIDE v2.5 22-Apr-2011 22:47:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @channelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @channelGUI_OutputFcn, ...
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


% --- Executes just before channelGUI is made visible.
function channelGUI_OpeningFcn(hObject, ~, handles, varargin)
%
% channelGUI('mainFig', handles.figure1) - call from movieDataGUI
% channelGUI(channelArray) - call from command line
% channelGUI(..., 'modal') - call channelGUI as a modal window
%
% User Data:
% 
% userData.channels - array of channel object
% userData.selectedChannel - index of current channel
% userData.mainFig - handle of movieDataGUI
% userData.properties - structure array of properties 
% userData.propNames - list of modifiable properties
%
% userData.helpFig - handle of help window
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for channelGUI
handles.output = hObject;

% Load help icon from dialogicons.mat
load lccbGuiIcons.mat
supermap(1,:) = get(hObject,'color');

userData.colormap = supermap;
userData.questIconData = questIconData;

axes(handles.axes_help);
Img = image(questIconData);
set(hObject,'colormap',supermap);
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
if openHelpFile
    set(Img, 'UserData', struct('class', mfilename))
end

if nargin > 3    
    if any(strcmp(varargin, 'mainFig'))
        % Called from movieDataGUI
        
        % Get userData.mainFig
        t = find(strcmp(varargin, 'mainFig'));
        userData.mainFig = varargin{t(end)+1};
        
        % Get userData.channels
        userData_main = get(userData.mainFig, 'UserData');
        assert(isa(userData_main.channels(1), 'Channel'), 'User-defined: No Channel object found.')
        userData.channels = userData_main.channels;
        
        % Get userData.selectedChannel
        handles_main = guidata(userData.mainFig);
        userData.selectedChannel = get(handles_main.listbox_channel, 'Value');
        
    elseif isa(varargin{1}(1), 'Channel')
        % Called from command line        
        userData.channels = varargin{1};
        userData.selectedChannel = 1;        
    else
        error('User-defined: Input parameters are incorrect.')
    end
    
    % Set as modal window
    if any(strcmp(varargin, 'modal'))
        set(hObject, 'WindowStyle', 'modal')
    end
    
else
    error('User-defined: No proper input.')
end

% Read channel initial properties and store them in a structure
userData.propNames={'excitationWavelength_','emissionWavelength_','exposureTime_'};
for i=1:numel(userData.channels)
    for j=1:numel(userData.propNames)
        propName= userData.propNames{j};
        userData.properties(i).(propName)=userData.channels(i).(propName);
    end
end

% Set up pop-up menu
set(handles.popupmenu_channel, 'String', ... 
arrayfun(@(x)(['Channel ' num2str(x)]), 1:length(userData.channels), 'UniformOutput', false) )
set(handles.popupmenu_channel, 'Value', userData.selectedChannel)

% Set up channel path and properties
set(handles.text_path, 'String', userData.channels(userData.selectedChannel).channelPath_);
for i=1:numel(userData.propNames)
    propHandle = handles.(['edit_' userData.propNames{i}(1:end-1)]);
    set(propHandle,'String',...
        num2str(userData.properties(userData.selectedChannel).(userData.propNames{i})));
end

% Update handles structure
set(handles.figure1,'UserData',userData)
guidata(hObject, handles);

% UIWAIT makes channelGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = channelGUI_OutputFcn(~, ~, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in popupmenu_channel.
function popupmenu_channel_Callback(hObject, ~, handles)

userData = get(handles.figure1, 'UserData');
if get(hObject,'Value') == userData.selectedChannel, return; end

% Retrieve handles and names of non-empty edit fields
propHandles = cellfun(@(x) handles.(['edit_' x(1:end-1)]),userData.propNames);
propStrings =get(propHandles,'String');
validProps = ~cellfun(@isempty,propStrings);
propNames=userData.propNames(validProps);
propHandles=propHandles(validProps);

% Save properties in the structure array
if ~isempty(propNames)
    for i=1:numel(propNames)
        userData.properties(userData.selectedChannel).(propNames{i})=...
            str2double(get(propHandles(i),'String'));
    end
end

% Update the selected channel path and properties
userData.selectedChannel=get(hObject,'Value'); 
set(handles.text_path, 'String', userData.channels(userData.selectedChannel).channelPath_)
for i=1:numel(userData.propNames)
    propHandle = handles.(['edit_' userData.propNames{i}(1:end-1)]);
    set(propHandle,'String',...
        num2str(userData.properties(userData.selectedChannel).(userData.propNames{i})));
end

set(handles.figure1,'UserData',userData)

% --- Executes when edit field is changed.
function edit_property_Callback(hObject, ~, handles)

if ~isempty(get(hObject,'String'))
    % Read property value
    propTag = get(hObject,'Tag');
    propName = [propTag(length('edit_')+1:end) '_'];
    propValue = str2double(get(hObject,'String'));

    % Test property value using the class static method
    if ~Channel.checkValue(propName,propValue)
        warndlg('Invalid property value','Setting Error','modal');
        set(hObject,'BackgroundColor',[1 .8 .8]);
        set(handles.popupmenu_channel, 'Enable','off');
        return
    end
end
set(hObject,'BackgroundColor',[1 1 1]);
set(handles.popupmenu_channel, 'Enable','on');

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)

delete(handles.figure1)

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

% Retrieve names and handles of non-empty fields
propHandles = cellfun(@(x) handles.(['edit_' x(1:end-1)]),userData.propNames);
propStrings =get(propHandles,'String');
validProps = ~cellfun(@isempty,propStrings);
propNames=userData.propNames(validProps);
propHandles=propHandles(validProps);

% Save properties in the structure array
if ~isempty(propNames)
    for i=1:numel(propNames)
        userData.properties(userData.selectedChannel).(propNames{i})=...
            str2double(get(propHandles(i),'String'));
    end
end

% Set properties to channel objects
for i=1:numel(userData.channels)
    set(userData.channels(i),userData.properties(i));
end

set(handles.figure1,'UserData',userData)
delete(handles.figure1)
        
% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end