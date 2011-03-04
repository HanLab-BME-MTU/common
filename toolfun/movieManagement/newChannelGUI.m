function varargout = newChannelGUI(varargin)
% NEWCHANNELGUI M-file for newChannelGUI.fig
%      NEWCHANNELGUI, by itself, creates a new NEWCHANNELGUI or raises the existing
%      singleton*.
%
%      H = NEWCHANNELGUI returns the handle to a new NEWCHANNELGUI or the handle to
%      the existing singleton*.
%
%      NEWCHANNELGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEWCHANNELGUI.M with the given input arguments.
%
%      NEWCHANNELGUI('Property','Value',...) creates a new NEWCHANNELGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before newChannelGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to newChannelGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help newChannelGUI

% Last Modified by GUIDE v2.5 04-Aug-2010 16:13:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @newChannelGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @newChannelGUI_OutputFcn, ...
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


% --- Executes just before newChannelGUI is made visible.
function newChannelGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%
% newChannelGUI('mainFig', handles.figure1) - call from newMovieDataGUI
% newChannelGUI(channelArray) - call from command line
% newChannelGUI(..., 'modal') - call newChannelGUI as a modal window
%
% User Data:
% 
% userData.channels - array of channel object
% userData.channelIndex - index of current channel
% userData.mainFig - handle of newMovieDataGUI
% userData.params - matrix of params 
% userData.numParams - the number of pamameters
%
% userData.helpFig - handle of help window
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for newChannelGUI
handles.output = hObject;

userData.numParams = 3;

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
    set(Img, 'UserData', struct('class', 'Channel'))
end

if nargin > 3
    
    if any(strcmp(varargin, 'mainFig'))
        % Called from newMovieDataGUI
        
        % Get userData.mainFig
        t = find(strcmp(varargin, 'mainFig'));
        userData.mainFig = varargin{t(end)+1};
        
        % Get userData.channels
        userData_main = get(userData.mainFig, 'UserData');
        assert(isa(userData_main.channels(1), 'Channel'), 'User-defined: No Channel object found.')
        userData.channels = userData_main.channels;
        
        % Get userData.channelIndex
        handles_main = guidata(userData.mainFig);
        userData.channelIndex = get(handles_main.listbox_channel, 'Value');
        
        % Get userData.params
        userData.params = userfcn_getParams(userData.channels, 1:length(userData.channels), userData.numParams);
    
    elseif isa(varargin{1}(1), 'Channel')
        % Called from command line
        
        userData.channels = varargin{1};
        userData.channelIndex = 1;
        userData.params = userfcn_getParams(userData.channels, 1:length(userData.channels), userData.numParams);
        
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

% Set up pop-up menu
set(handles.popupmenu_channel, 'String', ... 
arrayfun(@(x)(['Channel ' num2str(x)]), 1:length(userData.channels), 'UniformOutput', false) )
set(handles.popupmenu_channel, 'Value', userData.channelIndex(end))

% Set up channel path
set(handles.text_path, 'String', userData.channels(userData.channelIndex(end)).channelPath_);

% Set up channel parameters
set(handles.edit_1, 'String', num2str(userData.params{userData.channelIndex(end), 1}))
set(handles.edit_2, 'String', num2str(userData.params{userData.channelIndex(end), 2}))
set(handles.edit_3, 'String', num2str(userData.params{userData.channelIndex(end), 3}))


% Update handles structure
set(handles.figure1,'UserData',userData)
guidata(hObject, handles);

% UIWAIT makes newChannelGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = newChannelGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu_channel.
function popupmenu_channel_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

% Input validation
text1 = get(handles.edit_1, 'String');
text2 = get(handles.edit_2, 'String');
text3 = get(handles.edit_3, 'String');

if ~isempty(text1) && ( isnan(str2double(text1)) ...
   || str2double(text1) < 0 )
    errordlg('Please provide a valid input for ''Excitation Wavelength''.','Setting Error','modal');
    
    set(hObject, 'Value', userData.channelIndex(end))
    return
end
if ~isempty(text2) && ( isnan(str2double(text2)) ...
   || str2double(text2) < 0 )
    errordlg('Please provide a valid input for ''Emission Wavelength''.','Setting Error','modal');
    
    set(hObject, 'Value', userData.channelIndex(end))
    return
end
if ~isempty(text3) && ( isnan(str2double(text3)) ...
   || str2double(text3) < 0 )
    errordlg('Please provide a valid input for ''Exposure Time''.','Setting Error','modal');
    
    set(hObject, 'Value', userData.channelIndex(end))
    return
end

% Save parameters to userData.params
userData.params{userData.channelIndex(end), 1} = str2num(get(handles.edit_1, 'String'));
userData.params{userData.channelIndex(end), 2} = str2num(get(handles.edit_2, 'String'));
userData.params{userData.channelIndex(end), 3} = str2num(get(handles.edit_3, 'String'));

% Save channel index
userData.channelIndex = cat(2, userData.channelIndex, get(hObject, 'Value'));

% Set channel path and parameters
set(handles.text_path, 'String', userData.channels(userData.channelIndex(end)).channelPath_)
set(handles.edit_1, 'String', num2str(userData.params{userData.channelIndex(end), 1}))
set(handles.edit_2, 'String', num2str(userData.params{userData.channelIndex(end), 2}))
set(handles.edit_3, 'String', num2str(userData.params{userData.channelIndex(end), 3}))

set(handles.figure1,'UserData',userData)

% --- Executes during object creation, after setting all properties.
function popupmenu_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_3 as text
%        str2double(get(hObject,'String')) returns contents of edit_3 as a double


% --- Executes during object creation, after setting all properties.
function edit_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_4_Callback(hObject, eventdata, handles)
% hObject    handle to edit_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_4 as text
%        str2double(get(hObject,'String')) returns contents of edit_4 as a double


% --- Executes during object creation, after setting all properties.
function edit_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% 
userData = get(handles.figure1, 'UserData');

% Input validation
text1 = get(handles.edit_2, 'String');
text2 = get(handles.edit_2, 'String');
text3 = get(handles.edit_3, 'String');

if ~isempty(text1) && ( isnan(str2double(text1)) ...
   || str2double(text1) < 0 )
    errordlg('Please provide a valid input for ''Excitation Wavelength''.','Setting Error','modal');
    
    set(hObject, 'Value', userData.channelIndex(end))
    return
end
if ~isempty(text2) && ( isnan(str2double(text2)) ...
   || str2double(text2) < 0 )
    errordlg('Please provide a valid input for ''Emission Wavelength''.','Setting Error','modal');
    
    set(hObject, 'Value', userData.channelIndex(end))
    return
end
if ~isempty(text3) && ( isnan(str2double(text3)) ...
   || str2double(text3) < 0 )
    errordlg('Please provide a valid input for ''Exposure Time''.','Setting Error','modal');
    
    set(hObject, 'Value', userData.channelIndex(end))
    return
end

% Save parameters to userData.params
userData.params{userData.channelIndex(end), 1} = str2num(get(handles.edit_1, 'String'));
userData.params{userData.channelIndex(end), 2} = str2num(get(handles.edit_2, 'String'));
userData.params{userData.channelIndex(end), 3} = str2num(get(handles.edit_3, 'String'));

% Set params to channel objects
userData.channelIndex = unique(userData.channelIndex);
userfcn_setParams(userData.channels, userData.params, userData.channelIndex)

set(handles.figure1,'UserData',userData)
delete(handles.figure1)


function userfcn_setParams (channels, params, index)
%
l = length(channels);
arrayfun(@(x)assert( x <= l, 'User-defined: index is larger than the number of channels'), index)

for i = index 

    channels(i).setExcitationWavelength(params{i,1})
    channels(i).setEmissionWavelength(params{i,2})
    channels(i).setExposureTime(params{i,3})

end

function params = userfcn_getParams(channels, index, numParams)
%
l = length(channels);
arrayfun(@(x)assert( x <= l, 'User-defined: index is larger than the number of channels'), index)

params = cell(l, numParams);

for i = index
    
    params{i,1} = channels(i).excitationWavelength_;
    params{i,2} = channels(i).emissionWavelength_;
    params{i,3} = channels(i).exposureTime_;
    
end

        
    


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
