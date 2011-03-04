function varargout = thresholdProcessGUI(varargin)
% THRESHOLDPROCESSGUI M-file for thresholdProcessGUI.fig
%      THRESHOLDPROCESSGUI, by itself, creates a new THRESHOLDPROCESSGUI or raises the existing
%      singleton*.
%
%      H = THRESHOLDPROCESSGUI returns the handle to a new THRESHOLDPROCESSGUI or the handle to
%      the existing singleton*.
%
%      THRESHOLDPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHOLDPROCESSGUI.M with the given input arguments.
%
%      THRESHOLDPROCESSGUI('Property','Value',...) creates a new THRESHOLDPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before thresholdProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to thresholdProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help thresholdProcessGUI

% Last Modified by GUIDE v2.5 02-Sep-2010 14:49:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @thresholdProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @thresholdProcessGUI_OutputFcn, ...
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


% --- Executes just before thresholdProcessGUI is made visible.
function thresholdProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.set1Fig = thresholdProcessGUI('mainFig', handles.figure1, procID);
%
% Available tools 
% UserData data:
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main setting panel
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%       userData.procConstr - constructor of current process
%       userData.channelIndex - the array of index for channels
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for segmentationProcessGUI
handles.output = hObject;

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);
userData.channelIndex = get(userData.handles_main.listbox_2, 'UserData');

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.crtPackage = userData_main.crtPackage;

% Get current process constructer
userData.procConstr = userData_main.procConstr{userData.procID};

% Get current process
if ~isempty(userData_main.crtProc) && isa(userData_main.crtProc, userData_main.procName{userData.procID})
    userData.crtProc = userData_main.crtProc;
    
elseif ~isempty(userData_main.segProc{userData.procID})
    userData.crtProc = userData_main.segProc{userData.procID};
    
else
    % Create new process and handle the process to user data and
    % array of segmentation process in main setting panel    
    userData.crtProc = userData.procConstr(userData_main.MD, userData.crtPackage.outputDirectory_);
    userData_main.segProc{userData.procID} = userData.crtProc;
end


% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;


% ---------------------- Parameter Setup -------------------------

funParams = userData.crtProc.funParams_;

if isempty(funParams.ThresholdValue)
    if funParams.MaxJump
       set(handles.checkbox_max, 'Value', 1);
       set(handles.edit_jump, 'Enable', 'on', 'String',...
                num2str(funParams.MaxJump));
    end
else
    set(handles.text_body3, 'Enable', 'on')
    set(handles.edit_coef, 'Enable', 'on')
    set(handles.pushbutton_add, 'Enable', 'on')
    set(handles.pushbutton_up, 'Enable', 'on')
    set(handles.pushbutton_down, 'Enable', 'on')
    set(handles.pushbutton_coef_delete, 'Enable', 'on')
    set(handles.listbox_coef1, 'Enable', 'on')
    
    set(handles.checkbox_auto, 'Value', 0)
    set(handles.checkbox_max, 'Enable', 'off')
    
    
    threshold = cell(1, length(funParams.ThresholdValue));
    
    for i = 1:length(funParams.ThresholdValue)
       threshold{i} = funParams.ThresholdValue(i); 
    end
    
    set(handles.listbox_coef1, 'String', threshold)
end


% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
axes(handles.axes_help);
Img = image(userData.questIconData); 
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn);
if openHelpFile
    set(Img, 'UserData', struct('class',class(userData.crtProc)))
end

% ----------------------------------------------------------------

% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = thresholdProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');
delete(handles.figure1);

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button

userData = get(handles.figure1, 'UserData');

% -------- Check user input --------

if get(handles.checkbox_auto, 'value')
    if get(handles.checkbox_max, 'Value')
        % If both checkbox are checked
        if isnan(str2double(get(handles.edit_jump, 'String'))) ...
                || str2double(get(handles.edit_jump, 'String')) < 0
            errordlg('Please provide a valid input for ''Maximum threshold jump''.','Setting Error','modal');
            return;
        end    
    end
else
    threshold = get(handles.listbox_coef1, 'String');
    if isempty(threshold)
       errordlg('Please provide at least one threshold value.','Setting Error','modal')
       return
    elseif length(threshold) ~= 1 && length(threshold) ~= length(userData.channelIndex)
       errordlg('Please provide the same number of threshold values as the input channels.','Setting Error','modal')
       return
    else
        threshold = str2double(threshold);
        if any(isnan(threshold)) || any(threshold < 0)
            errordlg('Please provide valid threshold values. Threshold cannot be a negative number.','Setting Error','modal')
            return            
        end
    end
end
   
% -------- Set parameter --------

funParams = userData.crtProc.funParams_;
    
if userData.crtProc.procChanged_ 

    if get(handles.checkbox_auto, 'value')
        % if automatic thresholding
        funParams.ThresholdValue = [ ];
        if get(handles.checkbox_max, 'value')
            funParams.MaxJump = str2double(get(handles.edit_jump,'String'));
        else
            funParams.MaxJump = 0;
        end
    else
        % if fixed thresholding
        if length(threshold) == 1
            funParams.ThresholdValue = repmat(threshold, [1 length(userData.channelIndex)]);
        else
            funParams.ThresholdValue = threshold;
        end
    end
    % Set parameters
    userData.crtProc.setPara(funParams);
    
end

% Save user data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);




% --- Executes on button press in checkbox_auto.
function checkbox_auto_Callback(hObject, eventdata, handles)
switch get(hObject, 'Value')
    case 0
        set(handles.text_body3, 'Enable', 'on')
        set(handles.edit_coef, 'Enable', 'on')
        set(handles.pushbutton_add, 'Enable', 'on')
        set(handles.pushbutton_up, 'Enable', 'on')
        set(handles.pushbutton_down, 'Enable', 'on')
        set(handles.pushbutton_coef_delete, 'Enable', 'on')
        set(handles.listbox_coef1, 'Enable', 'on')
        
        
        set(handles.checkbox_max, 'Enable', 'off', 'Value', 0);
        set(handles.edit_jump, 'Enable', 'off');
       
    case 1
        set(handles.text_body3, 'Enable', 'off')
        set(handles.edit_coef, 'Enable', 'off')
        set(handles.pushbutton_add, 'Enable', 'off')
        set(handles.pushbutton_up, 'Enable', 'off')
        set(handles.pushbutton_down, 'Enable', 'off')
        set(handles.pushbutton_coef_delete, 'Enable', 'off')
        set(handles.listbox_coef1, 'Enable', 'off')        
        
        set(handles.checkbox_max, 'Enable', 'on');
                
        
end
userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);



function edit_jump_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes during object creation, after setting all properties.
function edit_jump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_max.
function checkbox_max_Callback(hObject, eventdata, handles)

switch get(hObject, 'value')
    case 0
        set(handles.edit_jump, 'Enable', 'off');
    case 1
        set(handles.edit_jump, 'Enable', 'on');
end

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes on selection change in listbox_coef1.
function listbox_coef1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_coef1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_coef1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_coef1


% --- Executes during object creation, after setting all properties.
function listbox_coef1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_coef1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_coef_delete.
function pushbutton_coef_delete_Callback(hObject, eventdata, handles)

% Call back function of 'delete' button
contents = get(handles.listbox_coef1,'String');
% Return if list is empty
if isempty(contents)
    return;
end
id = get(handles.listbox_coef1,'Value');

% Delete selected item
contents(id) = [ ];

% Refresh listbox
set(handles.listbox_coef1,'String',contents);
% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id>length(contents) && id>1)
    set(handles.listbox_coef1,'Value',length(contents));
end

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes on button press in pushbutton_up.
function pushbutton_up_Callback(hObject, eventdata, handles)
% % call back of 'Up' button
id = get(handles.listbox_coef1,'Value');
contents = get(handles.listbox_coef1,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == 1
    return;
end

temp = contents{id};
contents{id} = contents{id-1};
contents{id-1} = temp;

set(handles.listbox_coef1, 'string', contents);
set(handles.listbox_coef1, 'value', id-1);

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes on button press in pushbutton_down.
function pushbutton_down_Callback(hObject, eventdata, handles)
% Call back of 'Down' button
id = get(handles.listbox_coef1,'Value');
contents = get(handles.listbox_coef1,'String');

% Return if list is empty
if isempty(contents) || isempty(id) || id == length(contents)
    return;
end

temp = contents{id};
contents{id} = contents{id+1};
contents{id+1} = temp;

set(handles.listbox_coef1, 'string', contents);
set(handles.listbox_coef1, 'value', id+1);

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
set(handles.listbox_coef1, 'Value', 1);
text = get(handles.edit_coef, 'String');
if isempty(text)
    return;
end

if isnan(str2double(text)) || str2double(text) < 0 
    errordlg('Please provide a valid coefficient. Coefficient must be positive.','Setting Error','modal');
    return;
end

contents = get(handles.listbox_coef1, 'String');
contents{end + 1} = text;
set(handles.listbox_coef1, 'String', contents)
set(handles.edit_coef, 'String', '')

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);



function edit_coef_Callback(hObject, eventdata, handles)
% hObject    handle to edit_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_coef as text
%        str2double(get(hObject,'String')) returns contents of edit_coef as a double


% --- Executes during object creation, after setting all properties.
function edit_coef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_coef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
