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

% Last Modified by GUIDE v2.5 06-Mar-2013 17:15:12

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

set(handles.edit_subsample_number,'String',funParams.Sub_Sample_Num);

set(handles.edit_StPaceSize,'String',funParams.StPace_Size);
set(handles.edit_StPatchSize,'String',funParams.StPatch_Size);
set(handles.edit_st_lowerbound_localthresholding,'String',funParams.st_lowerbound_localthresholding);

set(handles.edit_IntPaceSize,'String',funParams.IntPace_Size);
set(handles.edit_IntPatchSize,'String',funParams.IntPatch_Size);
set(handles.edit_int_lowerbound_localthresholding,'String',funParams.int_lowerbound_localthresholding);

set(handles.popupmenu_cell_mask, 'Value',funParams.Cell_Mask_ind);

set(handles.checkbox_outgrowth,'value',funParams.VIF_Outgrowth_Flag);


Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only','geo_based'};
for Combine_Way_ind = 1 : 6
    if(strcmp(funParams.Combine_Way, Combine_Way_tag{Combine_Way_ind}))
        set(handles.popupmenu_segmentationbase,'Value',Combine_Way_ind);
        break;
    end
end


if (strcmp(funParams.Combine_Way,'st_only')||strcmp(funParams.Combine_Way,'st_nms_two')||strcmp(funParams.Combine_Way,'st_nms_only'))
    
    set(handles.edit_StPaceSize,'Enable','on');
    set(handles.edit_StPatchSize,'Enable','on');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
    set(handles.edit_IntPaceSize,'Enable','off');
    set(handles.edit_IntPatchSize,'Enable','off');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','off');
    
else
    if (strcmp(funParams.Combine_Way,'int_only'))
        
        set(handles.edit_StPaceSize,'Enable','off');
        set(handles.edit_StPatchSize,'Enable','off');
        set(handles.edit_st_lowerbound_localthresholding,'Enable','off');
        set(handles.edit_IntPaceSize,'Enable','on');
        set(handles.edit_IntPatchSize,'Enable','on');
        set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
        
    else
        set(handles.edit_StPaceSize,'Enable','on');
        set(handles.edit_StPatchSize,'Enable','on');
        set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
        set(handles.edit_IntPaceSize,'Enable','on');
        set(handles.edit_IntPatchSize,'Enable','on');
        set(handles.edit_int_lowerbound_localthresholding,'Enable','on');        
    end
end

if (strcmp(funParams.Combine_Way,'geo_based'))    
    set(handles.edit_StPaceSize,'Enable','off');
    set(handles.edit_StPatchSize,'Enable','off');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','off');
    set(handles.edit_IntPaceSize,'Enable','off');
    set(handles.edit_IntPatchSize,'Enable','off');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','off');
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

Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only','geo_based'};
Combine_Way_ind = get(handles.popupmenu_segmentationbase, 'Value');
funParams.Combine_Way=Combine_Way_tag{Combine_Way_ind};

StPace_Size = str2double(get(handles.edit_StPaceSize, 'String'));
if isnan(StPace_Size) || StPace_Size < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.StPace_Size=StPace_Size;

StPatch_Size = str2double(get(handles.edit_StPatchSize, 'String'));
if isnan(StPatch_Size) || StPatch_Size < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.StPatch_Size=StPatch_Size;


st_lowerbound_localthresholding = str2double(get(handles.edit_st_lowerbound_localthresholding, 'String'));
if isnan(st_lowerbound_localthresholding) || st_lowerbound_localthresholding < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.st_lowerbound_localthresholding=st_lowerbound_localthresholding;


IntPace_Size = str2double(get(handles.edit_IntPaceSize, 'String'));
if isnan(IntPace_Size) || IntPace_Size < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_PaceSize,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.IntPace_Size=IntPace_Size;

IntPatch_Size = str2double(get(handles.edit_IntPatchSize, 'String'));
if isnan(IntPatch_Size) || IntPatch_Size < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_Patch_Size,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.IntPatch_Size=IntPatch_Size;


int_lowerbound_localthresholding = str2double(get(handles.edit_int_lowerbound_localthresholding, 'String'));
if isnan(int_lowerbound_localthresholding) || int_lowerbound_localthresholding < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_lowerbound_localthresholding,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.int_lowerbound_localthresholding=int_lowerbound_localthresholding;


Cell_Mask_ind = get(handles.popupmenu_cell_mask, 'Value');
funParams.Cell_Mask_ind = Cell_Mask_ind;

VIF_Outgrowth_Flag = get(handles.checkbox_outgrowth,'value');
funParams.VIF_Outgrowth_Flag = VIF_Outgrowth_Flag;

Sub_Sample_Num  = str2double(get(handles.edit_subsample_number, 'String'));
if isnan(Sub_Sample_Num) || Sub_Sample_Num < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_subsample,'String') '''.'],'Setting Error','modal');
    return;
end

funParams.Sub_Sample_Num  = Sub_Sample_Num;
 

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

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_selectedChannels, 'String', contents2, 'Userdata', chanIndex2);


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

function edit_StPatchSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StPatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StPatchSize as text
%        str2double(get(hObject,'String')) returns contents of edit_StPatchSize as a double


% --- Executes during object creation, after setting all properties.
function edit_StPatchSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StPatchSize (see GCBO)
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



function edit_st_lowerbound_localthresholding_Callback(hObject, eventdata, handles)
% hObject    handle to edit_st_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_st_lowerbound_localthresholding as text
%        str2double(get(hObject,'String')) returns contents of edit_st_lowerbound_localthresholding as a double


% --- Executes during object creation, after setting all properties.
function edit_st_lowerbound_localthresholding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_st_lowerbound_localthresholding (see GCBO)
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


% The following code is for set the parameter set editing boxes, that corresponds to combine way 
% not choosen as disabled from editing

Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only','geo_based'};
Combine_Way_ind = get(hObject, 'Value');
Combine_Way=Combine_Way_tag{Combine_Way_ind};

if (strcmp(Combine_Way,'st_only')||strcmp(Combine_Way,'st_nms_two')||strcmp(Combine_Way,'st_nms_only'))
    set(handles.edit_StPaceSize,'Enable','on');
    set(handles.edit_StPatchSize,'Enable','on');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
    set(handles.edit_IntPaceSize,'Enable','off');
    set(handles.edit_IntPatchSize,'Enable','off');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','off');
    
else
    if (strcmp(Combine_Way,'int_only'))
        set(handles.edit_StPaceSize,'Enable','off');
        set(handles.edit_StPatchSize,'Enable','off');
        set(handles.edit_st_lowerbound_localthresholding,'Enable','off');
        set(handles.edit_IntPaceSize,'Enable','on');
        set(handles.edit_IntPatchSize,'Enable','on');
        set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
    else
        set(handles.edit_StPaceSize,'Enable','on');
        set(handles.edit_StPatchSize,'Enable','on');
        set(handles.edit_st_lowerbound_localthresholding,'Enable','on');
        set(handles.edit_IntPaceSize,'Enable','on');
        set(handles.edit_IntPatchSize,'Enable','on');
        set(handles.edit_int_lowerbound_localthresholding,'Enable','on');
    end
end


if (strcmp(Combine_Way,'geo_based'))    
    set(handles.edit_StPaceSize,'Enable','off');
    set(handles.edit_StPatchSize,'Enable','off');
    set(handles.edit_st_lowerbound_localthresholding,'Enable','off');
    set(handles.edit_IntPaceSize,'Enable','off');
    set(handles.edit_IntPatchSize,'Enable','off');
    set(handles.edit_int_lowerbound_localthresholding,'Enable','off');
end

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

set(hObject,'String',{'Steerable Filter Results','Intensity','Combine Both St and Int',...
    'St NMS two','St nms only', 'Geo based'});

Combine_Way_tag = {'st_only','int_only','int_st_both','st_nms_two','st_nms_only','geo_based'};
set(hObject,'Value',1);
% for Combine_Way_ind = 1 : 6
%     if(strcmp(funParams.Combine_Way, Combine_Way_tag{Combine_Way_ind}))
%         set(hObject,'Value',Combine_Way_ind);
%         break;
%     end
% end
 

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



function edit_subsample_number_Callback(hObject, eventdata, handles)
% hObject    handle to edit_subsample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_subsample_number as text
%        str2double(get(hObject,'String')) returns contents of edit_subsample_number as a double


% --- Executes during object creation, after setting all properties.
function edit_subsample_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_subsample_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IntPaceSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IntPaceSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IntPaceSize as text
%        str2double(get(hObject,'String')) returns contents of edit_IntPaceSize as a double


% --- Executes during object creation, after setting all properties.
function edit_IntPaceSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IntPaceSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_IntPatchSize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IntPatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IntPatchSize as text
%        str2double(get(hObject,'String')) returns contents of edit_IntPatchSize as a double


% --- Executes during object creation, after setting all properties.
function edit_IntPatchSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IntPatchSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_int_lowerbound_localthresholding_Callback(hObject, eventdata, handles)
% hObject    handle to edit_int_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_int_lowerbound_localthresholding as text
%        str2double(get(hObject,'String')) returns contents of edit_int_lowerbound_localthresholding as a double


% --- Executes during object creation, after setting all properties.
function edit_int_lowerbound_localthresholding_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_int_lowerbound_localthresholding (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_traing_nms.
function pushbutton_traing_nms_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_traing_nms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;
funParams.F_classifier = nms_classifier_train(userData.MD,[],[]);

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);


% --- Executes on button press in pushbutton_load_nms_classifier.
function pushbutton_load_nms_classifier_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_nms_classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

funParams.F_classifier = load_classifier_trained(userData.MD);


% --- Executes during object creation, after setting all properties.
function pushbutton_traing_nms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_traing_nms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
