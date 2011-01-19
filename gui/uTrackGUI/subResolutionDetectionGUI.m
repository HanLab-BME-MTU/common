function varargout = subResolutionDetectionGUI(varargin)
% SUBRESOLUTIONDETECTIONGUI M-file for subResolutionDetectionGUI.fig
%      SUBRESOLUTIONDETECTIONGUI, by itself, creates a new SUBRESOLUTIONDETECTIONGUI or raises the existing
%      singleton*.
%
%      H = SUBRESOLUTIONDETECTIONGUI returns the handle to a new SUBRESOLUTIONDETECTIONGUI or the handle to
%      the existing singleton*.
%
%      SUBRESOLUTIONDETECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SUBRESOLUTIONDETECTIONGUI.M with the given input arguments.
%
%      SUBRESOLUTIONDETECTIONGUI('Property','Value',...) creates a new SUBRESOLUTIONDETECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before subResolutionDetectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to subResolutionDetectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help subResolutionDetectionGUI

% Last Modified by GUIDE v2.5 08-Dec-2010 17:49:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @subResolutionDetectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @subResolutionDetectionGUI_OutputFcn, ...
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


% --- Executes just before subResolutionDetectionGUI is made visible.
function subResolutionDetectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.set1Fig = subResolutionDetectionGUI('mainFig', handles.figure1, procID);
%

% Available tools 
% UserData data:
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main setting panel
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%       userData.procConstr - constructor of current process
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%
%       userData.path
%       userData.file          
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

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.MD = userData_main.MD;
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

% handles update
handles.edit_av = [handles.edit_av_1 handles.edit_av_2 handles.edit_av_3];
handles.edit_wz = [handles.edit_ws_1 handles.edit_ws_2 handles.edit_ws_3];

% ---------------------- Parameter Setup -------------------------

funParams = userData.crtProc.funParams_;

% funParams.detectionParam

set(handles.edit_gsd, 'String', num2str(funParams.detectionParam.psfSigma))
set(handles.edit_cbd, 'String', num2str(funParams.detectionParam.bitDepth))

arrayfun(@(x)set(handles.edit_av(x), 'String', num2str(funParams.detectionParam.alphaLocMax(x))), ...
                  1:length(funParams.detectionParam.alphaLocMax))
              
if  length(funParams.detectionParam.integWindow)>1 || funParams.detectionParam.integWindow
    
    set(handles.checkbox_rollingwindow, 'Value', 1)
    set(handles.text_body_6, 'Enable', 'on')
    arrayfun(@(x)set(handles.edit_wz(x), 'Enable', 'on'), 1:3)
    arrayfun(@(x)set(handles.edit_wz(x), 'String', num2str(funParams.detectionParam.integWindow(x)*2+1)), ...
                  1:length(funParams.detectionParam.integWindow))
end

if funParams.detectionParam.doMMF
    
    set(handles.checkbox_mmf, 'Value', 1)
    set(handles.text_body_8, 'Enable', 'on')
    set(handles.edit_r, 'Enable', 'on')
end

set(handles.edit_r, 'String', num2str(funParams.detectionParam.testAlpha.alphaR))
set(handles.edit_a, 'String', num2str(funParams.detectionParam.testAlpha.alphaA))
set(handles.edit_d, 'String', num2str(funParams.detectionParam.testAlpha.alphaD))
set(handles.edit_f, 'String', num2str(funParams.detectionParam.testAlpha.alphaF))

if funParams.detectionParam.numSigmaIter
   
    set(handles.checkbox_iteration, 'Value', 1)
    set(handles.text_body_12, 'Enable', 'on')
    set(handles.edit_num, 'Enable', 'on', 'String',num2str(funParams.detectionParam.numSigmaIter) )
end

if funParams.detectionParam.visual
    
   set(handles.checkbox_visual, 'Value', 1) 
end



if ~isempty(funParams.detectionParam.background)
    
    set(handles.checkbox_background, 'Value', 1)
    checkbox_background_Callback(handles.checkbox_background, [], handles)
    
    set(handles.edit_av_background, 'String', num2str(funParams.detectionParam.background.alphaLocMaxAbs))
    set(handles.edit_path, 'String', funParams.detectionParam.background.imageDir)

end

% funParams.movieParam

set(handles.edit_min, 'String', num2str(funParams.movieParam.firstImageNum))
set(handles.edit_max, 'String', num2str(funParams.movieParam.lastImageNum))
set(handles.text_framenum, 'String', ['(Totally ' num2str(userData.crtProc.owner_.nFrames_) ' frames in the movie)'])

% funParams.saveResults

userData.path = funParams.saveResults.dir;
userData.file = userData.crtProc.filename_;

str = [funParams.saveResults.dir userData.crtProc.filename_];
if length(str)>100
    str = ['... ' str(end-500:end)];
end
set(handles.text_path, 'String', str)


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
    set(Img, 'UserData', struct('class', 'SubResolutionProcess'))
else
    set(Img, 'UserData', userData.crtProc.getHelp(true))
end

% ----------------------------------------------------------------

% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = subResolutionDetectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_path_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_path_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_path_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_saveas.
function pushbutton_saveas_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

% Ask user where to save the movie data file
[file path] = uiputfile('*.mat','Save Results As',...
             [userData.path userData.file]);
        
if ~any([file path])
    return;
end

userData.path = path;
userData.file = file;

str = [path file];
if length(str)>100
    str = ['... ' str(end-50:end)];
end
set(handles.text_path, 'String', str)

set(handles.figure1, 'UserData', userData)




function edit_min_Callback(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_min as text
%        str2double(get(hObject,'String')) returns contents of edit_min as a double


% --- Executes during object creation, after setting all properties.
function edit_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_max_Callback(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_max as text
%        str2double(get(hObject,'String')) returns contents of edit_max as a double


% --- Executes during object creation, after setting all properties.
function edit_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gsd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gsd as text
%        str2double(get(hObject,'String')) returns contents of edit_gsd as a double


% --- Executes during object creation, after setting all properties.
function edit_gsd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gsd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_cbd_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cbd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cbd as text
%        str2double(get(hObject,'String')) returns contents of edit_cbd as a double


% --- Executes during object creation, after setting all properties.
function edit_cbd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cbd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_mmf.
function checkbox_mmf_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_mmf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_mmf

if get(hObject,'Value')
   
    set(handles.text_body_8, 'Enable', 'on')
    set(handles.edit_r, 'Enable', 'on')
else
    
    set(handles.text_body_8, 'Enable', 'off')
    set(handles.edit_r, 'Enable', 'off')
end


function edit_r_Callback(hObject, eventdata, handles)
% hObject    handle to edit_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_r as text
%        str2double(get(hObject,'String')) returns contents of edit_r as a double


% --- Executes during object creation, after setting all properties.
function edit_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_f_Callback(hObject, eventdata, handles)
% hObject    handle to edit_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_f as text
%        str2double(get(hObject,'String')) returns contents of edit_f as a double


% --- Executes during object creation, after setting all properties.
function edit_f_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_d_Callback(hObject, eventdata, handles)
% hObject    handle to edit_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_d as text
%        str2double(get(hObject,'String')) returns contents of edit_d as a double


% --- Executes during object creation, after setting all properties.
function edit_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_a_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_a as text
%        str2double(get(hObject,'String')) returns contents of edit_a as a double


% --- Executes during object creation, after setting all properties.
function edit_a_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_iteration.
function checkbox_iteration_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_iteration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_iteration
if get(hObject, 'Value')
   
    set(handles.text_body_12, 'Enable', 'on')
    set(handles.edit_num, 'Enable', 'on')
    
    % If no input, use default number 10 
    if isempty(get(handles.edit_num, 'String')) ||...
            isnan(str2double(get(handles.edit_num, 'String')))
       
        set(handles.edit_num, 'String', '10')
    end

else
   
    set(handles.text_body_12, 'Enable', 'off')
    set(handles.edit_num, 'Enable', 'off')

    
end


function edit_num_Callback(hObject, eventdata, handles)
% hObject    handle to edit_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_num as text
%        str2double(get(hObject,'String')) returns contents of edit_num as a double


% --- Executes during object creation, after setting all properties.
function edit_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_rollingwindow.
function checkbox_rollingwindow_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_rollingwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_rollingwindow

if get(hObject, 'Value')
   
    set(handles.text_body_6, 'Enable', 'on')
    set(handles.edit_ws_1, 'Enable', 'on')
    set(handles.edit_ws_2, 'Enable', 'on')
    set(handles.edit_ws_3, 'Enable', 'on')
    
    % If no input, use default number 1
    if ~isempty(get(handles.edit_av_1, 'String')) && ...
            (isempty(get(handles.edit_ws_1, 'String')) ||...
            isnan(str2double(get(handles.edit_ws_1, 'String'))) )
       
        set(handles.edit_ws_1, 'String', '1')
    end
    
    if ~isempty(get(handles.edit_av_2, 'String')) && ...
            (isempty(get(handles.edit_ws_2, 'String')) ||...
            isnan(str2double(get(handles.edit_ws_2, 'String'))) )
       
        set(handles.edit_ws_2, 'String', '1')
    end
    
    if ~isempty(get(handles.edit_av_3, 'String')) && ...
            (isempty(get(handles.edit_ws_3, 'String')) ||...
            isnan(str2double(get(handles.edit_ws_3, 'String'))) )
       
        set(handles.edit_ws_3, 'String', '1')
    end    
    
else
   
    set(handles.text_body_6, 'Enable', 'off')
    set(handles.edit_ws_1, 'Enable', 'off')
    set(handles.edit_ws_2, 'Enable', 'off')
    set(handles.edit_ws_3, 'Enable', 'off')    
    
end



function edit_ws_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ws_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ws_1 as text
%        str2double(get(hObject,'String')) returns contents of edit_ws_1 as a double


% --- Executes during object creation, after setting all properties.
function edit_ws_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ws_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_path_Callback(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_path as text
%        str2double(get(hObject,'String')) returns contents of edit_path as
%        a double


% --- Executes during object creation, after setting all properties.
function edit_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_new.
function pushbutton_new_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_av_background_Callback(hObject, eventdata, handles)
% hObject    handle to edit_av_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_av_background as text
%        str2double(get(hObject,'String')) returns contents of edit_av_background as a double


% --- Executes during object creation, after setting all properties.
function edit_av_background_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_av_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_open.
function pushbutton_open_Callback(hObject, eventdata, handles)

pathname = uigetdir(pwd);
if isnumeric(pathname)
    return;
end

set(handles.edit_path, 'String', pathname)





% --- Executes on button press in checkbox_background.
function checkbox_background_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_background
if get(hObject, 'Value')
   
    set(handles.text_body_13, 'Enable', 'on')
    set(handles.edit_av_background, 'Enable', 'on')
    set(handles.edit_path, 'Enable', 'on')
    set(handles.pushbutton_new, 'Enable', 'on')
    set(handles.pushbutton_open, 'Enable', 'on')
    
else
   
    set(handles.text_body_13, 'Enable', 'off')
    set(handles.edit_av_background, 'Enable', 'off')
    set(handles.edit_path, 'Enable', 'off')
    set(handles.pushbutton_new, 'Enable', 'off')
    set(handles.pushbutton_open, 'Enable', 'off')  
    
end

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button

userData = get(handles.figure1, 'UserData');

% ------------------------- Check user input -----------------------------

missPara = [];

%% Required Fields

psfSigma = ( get(handles.edit_gsd, 'String') );
bitDepth = ( get(handles.edit_cbd, 'String') );

alphaLocMax{1} = ( get(handles.edit_av_1, 'String') );
alphaLocMax{2} = ( get(handles.edit_av_2, 'String') );
alphaLocMax{3} = ( get(handles.edit_av_3, 'String') );
integWindow{1} = ( get(handles.edit_ws_1, 'String') );
integWindow{2} = ( get(handles.edit_ws_2, 'String') );
integWindow{3} = ( get(handles.edit_ws_3, 'String') );

alphaR = ( get(handles.edit_r, 'String') );
alphaA = ( get(handles.edit_a, 'String') );
alphaD = ( get(handles.edit_d, 'String') );
alphaF = ( get(handles.edit_f, 'String') );

numSigmaIter = ( get(handles.edit_num, 'String') );
alphaLocMaxAbs = ( get(handles.edit_av_background, 'String') );
bg_dir = get(handles.edit_path, 'String');

minid = ( get(handles.edit_min, 'String') );
maxid = ( get(handles.edit_max, 'String') );

%% Gaussian Standard Deviation: psfSigma
if isempty( psfSigma )
    errordlg('Parameter "Gaussian Standard Deviation" is requied by the algorithm.','Error','modal')
    return
    
elseif isnan(str2double(psfSigma)) || str2double(psfSigma) < 0
    errordlg('Please provide a valid value to parameter "Gaussian Standard Deviation".','Error','modal')
    return
else
    psfSigma = str2double(psfSigma);
end

%% Camera Bit Depth: bitDepth
if isempty( bitDepth )
    errordlg('Parameter "Camera Bit Depth" is requied by the algorithm.','Error','modal')
    return
    
elseif isnan(str2double(bitDepth)) || str2double(bitDepth) < 0
    errordlg('Please provide a valid value to parameter "Camera Bit Depth".','Error','modal')
    return
    
else
    bitDepth = str2double(bitDepth);
end

%% Alpha-value for Local Maxima Detection: alphaLoc
temp = cellfun(@(x)~isempty(x), alphaLocMax);
if all(~temp)
    missPara = strcat(missPara, sprintf('Alpha-value for Local Maxima Detection\n'));
    alphaLoc = [];
    
elseif any(cellfun(@(x)isnan(str2double(x)), alphaLocMax(temp))) || any(cellfun(@(x)(str2double(x) < 0), alphaLocMax(temp)))
    errordlg('Please provide a valid value to parameter "Alpha-value for Local Maxima Detection".','Error','modal')
    return
    
else
    alphaLoc = cellfun(@(x)str2double(x), alphaLocMax(temp), 'UniformOutput', true);
end

%% Window Size: integWin
temp = cellfun(@(x)~isempty(x), integWindow);
if get(handles.checkbox_rollingwindow, 'Value')
    if isempty(alphaLoc)
        missPara = strcat(missPara, sprintf('Window Size\n'));
        integWin = (1-1)/2;
        
    elseif length(find(temp)) ~= length(alphaLoc)
        errordlg('The length of parameter "Camera Bit Depth" must be the same with parameter "Alpha-value for Local Maxima Detection".','Error','modal')
        return
        
    elseif any(cellfun(@(x)isnan(str2double(x)), integWindow(temp))) || any(cellfun(@(x)(str2double(x) < 0), integWindow(temp)))
        errordlg('Please provide a valid value to parameter "Camera Bit Depth".','Error','modal')   
        return
        
    else
        integWin = cellfun(@(x)str2double(x), integWindow(temp), 'UniformOutput', true);
        if ~all(mod(integWin, 2) == 1)
            errordlg('Parameter "Window Size" must be an odd number.','Error','modal')   
            return            
        end
        integWin = (integWin-1)/2;
    end
    
else
    integWin = 0;
end


%% AlphaR

if get(handles.checkbox_mmf, 'Value')
    if isempty( alphaR )
        missPara = strcat(missPara, sprintf('Alpha Residual\n'));
        alphaR = []; 
        
    elseif isnan(str2double(alphaR)) || str2double(alphaR) < 0
        errordlg('Please provide a valid value to parameter "Alpha Residuals".','Error','modal')
        return   
        
    else
        alphaR = str2double(alphaR);
    end
else
    alphaR = [];
end

%% AlphaA

if isempty( alphaA )
    missPara = strcat(missPara, sprintf('Alpha Amplitude\n'));
    alphaA = []; 
        
elseif isnan(str2double(alphaA)) || str2double(alphaA) < 0
    errordlg('Please provide a valid value to parameter "Alpha Amplitude".','Error','modal')
    return   
        
else
    alphaA = str2double(alphaA);
end

%% AlphaD

if isempty( alphaD )
    missPara = strcat(missPara, sprintf('Alpha Distance\n'));
    alphaD = []; 
        
elseif isnan(str2double(alphaD)) || str2double(alphaD) < 0
    errordlg('Please provide a valid value to parameter "Alpha Distance".','Error','modal')
    return   
        
else
    alphaD = str2double(alphaD);
end

%% AlphaF

if isempty( alphaF )
    missPara = strcat(missPara, sprintf('Alpha Final\n'));
    alphaD = []; 
        
elseif isnan(str2double(alphaF)) || str2double(alphaF) < 0
    errordlg('Please provide a valid value to parameter "Alpha Final".','Error','modal')
    return   
        
else
    alphaF = str2double(alphaF);
end

%% Maximum Number of Interations: numSigmaIter

if get(handles.checkbox_iteration, 'Value')
    if isempty( numSigmaIter )
        missPara = strcat(missPara, sprintf('Maximum Number of Interations\n'));
        numSigmaIter = []; 

    elseif isnan(str2double(numSigmaIter)) || str2double(numSigmaIter) < 0 ...
            || floor(str2double(numSigmaIter)) ~= ceil(str2double(numSigmaIter))
        errordlg('Please provide a valid value to parameter "Maximum Number of Interations".','Error','modal')
        return   

    else
        numSigmaIter = str2double(numSigmaIter);
    end    
    
else
    numSigmaIter = 0;
end

%% Alpha-value for Local Maxima Detection

if get(handles.checkbox_background, 'Value')
    if isempty( alphaLocMaxAbs )
        missPara = strcat(missPara, sprintf('Background Alpha-value \n'));
        alphaLocMaxAbs = []; 

    elseif isnan(str2double(alphaLocMaxAbs)) || str2double(alphaLocMaxAbs) < 0 
        errordlg('Please provide a valid value to parameter "Background Alpha-value".','Error','modal')
        return   

    else
        alphaLocMaxAbs = str2double(alphaLocMaxAbs);
    end    
    
else
    alphaLocMaxAbs = [];
end

%% Background directory: bg_dir

if get(handles.checkbox_background, 'Value')
    if isempty( bg_dir )
        errordlg('Please specify a background image directory.','Error','modal')
        return
    end    
    
else
    bg_dir = [];
end

%% Frame Index: minid, maxid

if isempty( minid )
    errordlg('Please specify the minimum index of frame.','Error','modal')
    return
        
elseif isnan(str2double(minid)) || str2double(minid) < 0 ...
        || floor(str2double(minid)) ~= ceil(str2double(minid))
    errordlg('Please provide a valid value to parameter "Minimum Frame Index".','Error','modal')
    return   
        
else
    minid = str2double(minid);
end

if isempty( maxid )
    errordlg('Please specify the maximum index of frame.','Error','modal')
    return
        
elseif isnan(str2double(maxid)) || str2double(maxid) < 0 ...
        || floor(str2double(maxid)) ~= ceil(str2double(maxid))
    errordlg('Please provide a valid value to parameter "Maximum Frame Index".','Error','modal')
    return   
        
else
    maxid = str2double(maxid);
end

if minid > maxid
    errordlg('Minimum frame index is larger than maximum frame index.','Error','modal')
    return       
end


% psfSigma, bitDepth, alphaLoc,integWin,alphaR,alphaA,alphaD,alphaF,numSigmaIter,alphaLocMaxAbs,bg_dir,minid,maxid

%% TO-DO check validation of output dir and background image dir

% Check background image directory
if get(handles.checkbox_background, 'Value')
    
    fileNames = imDir(bg_dir);
    if isempty(fileNames)
        errordlg(sprintf('No valid image file found in the background image directory: %s.', bg_dir),'Error','modal')
        return           
    end

    [x1 filenameBase x3 x4] =  getFilenameBody(fileNames(1).name);
end

%% ------------------------- Set Parameter -------------------------------

funParams = userData.crtProc.funParams_;

%% funParams.movieParam

funParams.movieParam.firstImageNum = minid;
funParams.movieParam.lastImageNum = maxid;

%% funParams.detectionParam

funParams.detectionParam.psfSigma = psfSigma;
funParams.detectionParam.bitDepth = bitDepth;
funParams.detectionParam.alphaLocMax = alphaLoc;
funParams.detectionParam.integWindow = integWin;

if get(handles.checkbox_mmf, 'Value')
    funParams.detectionParam.doMMF = 1;
else
    funParams.detectionParam.doMMF = 0;
end

funParams.detectionParam.testAlpha.alphaR = alphaR;
funParams.detectionParam.testAlpha.alphaA = alphaA;
funParams.detectionParam.testAlpha.alphaD = alphaD;
funParams.detectionParam.testAlpha.alphaF = alphaF;

funParams.detectionParam.numSigmaIter = numSigmaIter;

if get(handles.checkbox_visual, 'Value')
    funParams.detectionParam.visual = 1;
else
    funParams.detectionParam.visual = 0;
end

if get(handles.checkbox_background, 'Value')
    
    % If background image dir is changed, set process "changed"
    if ~(isfield(funParams.detectionParam.background, 'imageDir') && ...
            strcmp(funParams.detectionParam.background.imageDir, bg_dir))

        userData.crtProc.setProcChanged(true);
    end
    
    funParams.detectionParam.background.imageDir = bg_dir;
    funParams.detectionParam.background.alphaLocMaxAbs = alphaLocMaxAbs;
    funParams.detectionParam.background.filenameBase = filenameBase;
else
    funParams.detectionParam.background = [];
end

%% funParams.saveResults

funParams.saveResults.dir = userData.path;

%% Save Parameters

% Save funParams
userData.crtProc.setPara(funParams);

% Save result file name
userData.crtProc.setFileName(userData.file)


% Save user data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% Notice user if any parameters are missing
if ~isempty(missPara)
    
    text = sprintf('The following pamameters are not given and will use default values.\n\nDefaults:\n%s \n\nDo you want to finish setting and use the dafault values?', missPara);
    user_response = questdlg(text, 'Default Settings', 'Cancel', 'Finish', 'Finish');
    
    if strcmpi(user_response, 'cancel')
        return
    end    
end

delete(handles.figure1);


function edit_av_1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_av_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_av_1 as text
%        str2double(get(hObject,'String')) returns contents of edit_av_1 as a double


% --- Executes during object creation, after setting all properties.
function edit_av_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_av_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_av_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_av_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_av_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_av_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_av_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_av_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_av_3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_av_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_av_3 as text
%        str2double(get(hObject,'String')) returns contents of edit_av_3 as a double


% --- Executes during object creation, after setting all properties.
function edit_av_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_av_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ws_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ws_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ws_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_ws_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_ws_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ws_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ws_3_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ws_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ws_3 as text
%        str2double(get(hObject,'String')) returns contents of edit_ws_3 as a double


% --- Executes during object creation, after setting all properties.
function edit_ws_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ws_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_visual.
function checkbox_visual_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_visual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_visual


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
