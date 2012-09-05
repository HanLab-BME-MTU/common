function varargout = filamentSegmentationProcessGUI(varargin)
% filamentSegmentationProcessGUI M-file for filamentSegmentationProcessGUI.fig
%      filamentSegmentationProcessGUI, by itself, creates a new filamentSegmentationProcessGUI or raises the existing
%      singleton*.
%
%      H = filamentSegmentationProcessGUI returns the handle to a new filamentSegmentationProcessGUI or the handle to
%      the existing singleton*.
%
%      filamentSegmentationProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in filamentSegmentationProcessGUI.M with the given input arguments.
%
%      filamentSegmentationProcessGUI('Property','Value',...) creates a new filamentSegmentationProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before filamentSegmentationProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to filamentSegmentationProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help filamentSegmentationProcessGUI

% Last Modified by GUIDE v2.5 13-Aug-2012 15:59:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @filamentSegmentationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @filamentSegmentationProcessGUI_OutputFcn, ...
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


% --- Executes just before filamentSegmentationProcessGUI is made visible.
function filamentSegmentationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',0);

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableChannels,'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));

channelIndex = funParams.ChannelIndex;

% Find any parent process
userData.parentProc = userData.crtPackage.getParent(userData.procID);
if isempty(userData.crtPackage.processes_{userData.procID}) && ~isempty(userData.parentProc)
    % Check existence of all parent processes
    emptyParentProc = any(cellfun(@isempty,userData.crtPackage.processes_(userData.parentProc)));
    if ~emptyParentProc
        % Intersect channel index with channel index of parent processes
        parentChannelIndex = @(x) userData.crtPackage.processes_{x}.funParams_.ChannelIndex;
        for i = userData.parentProc
            channelIndex = intersect(channelIndex,parentChannelIndex(i));
        end
    end
   
end

if ~isempty(channelIndex)
    channelString = userData.MD.getChannelPaths(channelIndex);
else
    channelString = {};
end

set(handles.listbox_selectedChannels,'String',channelString,...
    'UserData',channelIndex);

set(handles.edit_PaceSize,'String',funParams.Pace_Size);
set(handles.edit_PatchSize,'String',funParams.Patch_Size);

set(handles.edit_lowerbound_localthresholding,'String',funParams.lowerbound_localthresholding);

set(handles.popupmenu_cell_mask, 'Value',funParams.Cell_Mask_ind);

set(handles.checkbox_outgrowth,'value',funParams.VIF_Outgrowth_Flag);

if (strcmp(funParams.Combine_Way,'st_only'))
    set(handles.popupmenu_segmentationbase, 'Value',1);
else
    if (strcmp(funParams.Combine_Way,'int_only'))
        set(handles.popupmenu_segmentationbase, 'Value',2);
    else
        set(handles.popupmenu_segmentationbase, 'Value',3);
    end
end
    
% Update user data and GUI data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = filamentSegmentationProcessGUI_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
end
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

Combine_Way_tag = {'st_only','int_only','int_st_both'};
Combine_Way_ind = get(handles.popupmenu_segmentationbase, 'Value');
funParams.Combine_Way=Combine_Way_tag{Combine_Way_ind};

Pace_Size = str2double(get(handles.edit_PaceSize, 'String'));
if isnan(Pace_Size) || Pace_Size < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.Pace_Size=Pace_Size;

Patch_Size = str2double(get(handles.edit_PatchSize, 'String'));
if isnan(Patch_Size) || Patch_Size < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.Patch_Size=Patch_Size;


lowerbound_localthresholding = str2double(get(handles.edit_lowerbound_localthresholding, 'String'));
if isnan(lowerbound_localthresholding) || lowerbound_localthresholding < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.lowerbound_localthresholding=lowerbound_localthresholding;

Cell_Mask_ind = get(handles.popupmenu_cell_mask, 'Value');
funParams.Cell_Mask_ind = Cell_Mask_ind;

VIF_Outgrowth_Flag = get(handles.checkbox_outgrowth,'value');
funParams.VIF_Outgrowth_Flag = VIF_Outgrowth_Flag;

funParams.OutputDirectory  = [ userData.crtPackage.outputDirectory_, filesep 'FilamentSegmentation'];

for iChannel = channelIndex
FilamentSegmentationOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(FilamentSegmentationOutputDir,'dir'))
        mkdir(FilamentSegmentationOutputDir);
    end
    
    userData.crtProc.setOutImagePath(iChannel,FilamentSegmentationOutputDir)
end   

% -------- Process Sanity check --------
% ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
contents1 = get(handles.listbox_availableChannels, 'String');

chanIndex1 = get(handles.listbox_availableChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedChannels, 'String', contents1);
        chanIndex2 = chanIndex1;
        thresholdValues =zeros(1,numel(chanIndex1));
    case 0
        set(handles.listbox_selectedChannels, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
        thresholdValues = [];
end
set(handles.listbox_selectedChannels, 'UserData', chanIndex2);
set(handles.listbox_thresholdValues,'String',num2cell(thresholdValues))

% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button

contents1 = get(handles.listbox_availableChannels, 'String');
contents2 = get(handles.listbox_selectedChannels, 'String');
id = get(handles.listbox_availableChannels, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_availableChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');
thresholdValues = cellfun(@str2num,get(handles.listbox_thresholdValues,'String'));

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        thresholdValues(end+1) = 0;
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_selectedChannels, 'String', contents2, 'Userdata', chanIndex2);
set(handles.listbox_thresholdValues,'String',num2cell(thresholdValues))


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
contents = get(handles.listbox_selectedChannels,'String');
id = get(handles.listbox_selectedChannels,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_selectedChannels, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_selectedChannels,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_selectedChannels,'String',contents);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if ishandle(userData.helpFig), delete(userData.helpFig); end
if ishandle(userData.previewFig), delete(userData.previewFig); end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

function edit_PatchSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PatchSize as text
%        str2double(get(hObject,'String')) returns contents of edit_PatchSize as a double


% --- Executes during object creation, after setting all properties.
function edit_PatchSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_lowerbound_localthresholding_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lowerbound_localthresholding as text
%        str2double(get(hObject,'String')) returns contents of edit_lowerbound_localthresholding as a double


% --- Executes during object creation, after setting all properties.
function edit_lowerbound_localthresholding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_segmentationbase.
function popupmenu_segmentationbase_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_segmentationbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_segmentationbase contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_segmentationbase


% --- Executes during object creation, after setting all properties.
function popupmenu_segmentationbase_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_segmentationbase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',{'Steerable Filter Results','Intensity','Combine Both'});
set(hObject,'Value',1);
 

% --- Executes on button press in checkbox_outgrowth.
function checkbox_outgrowth_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_outgrowth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_outgrowth


% --- Executes on selection change in popupmenu_cell_mask.
function popupmenu_cell_mask_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_cell_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_cell_mask contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_cell_mask


% --- Executes during object creation, after setting all properties.
function popupmenu_cell_mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_cell_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',{'Cell Segmentation in this same channel','Input ROI','Cell Segmentation combined from two channels','No limitation'});
set(hObject,'Value',4);
