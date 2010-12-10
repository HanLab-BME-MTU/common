function varargout = maskRefineProcessGUI(varargin)
% MASKREFINEPROCESSGUI M-file for maskRefineProcessGUI.fig
%      MASKREFINEPROCESSGUI, by itself, creates a new MASKREFINEPROCESSGUI or raises the existing
%      singleton*.
%
%      H = MASKREFINEPROCESSGUI returns the handle to a new MASKREFINEPROCESSGUI or the handle to
%      the existing singleton*.
%
%      MASKREFINEPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASKREFINEPROCESSGUI.M with the given input arguments.
%
%      MASKREFINEPROCESSGUI('Property','Value',...) creates a new MASKREFINEPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before maskRefineProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to maskRefineProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help maskRefineProcessGUI

% Last Modified by GUIDE v2.5 08-Sep-2010 15:28:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @maskRefineProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @maskRefineProcessGUI_OutputFcn, ...
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


% --- Executes just before maskRefineProcessGUI is made visible.
function maskRefineProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.set1Fig = maskRefineProcessGUI('mainFig', handles.figure1);
%
% Available tools 
% UserData data:
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main setting panel
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%       userData.procConstr - constructor of current process
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%
%       userData.helpFig - help figure
%

set(handles.text_copyright, 'String', userfcn_copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for segmentationProcessGUI
handles.output = hObject;

% Get main figure handle
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.handles_main = guidata(userData.mainFig);

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.crtPackage = userData_main.crtPackage;

% Get current process constructer
userData.procConstr = @MaskRefinementProcess;

% Get current process
if isempty(userData_main.maskRefineProc)
    % If no mask refinement process on GUI
    
    if ~isempty(userData_main.crtProc) && ~isempty(userData_main.crtProc.maskRefineProcess_)...
            && isa(userData_main.crtProc.maskRefineProcess_, 'MaskRefinementProcess')
        % If segmentation process has a mask refinement process attached to it,
        % then create a new mask refinement process with the same parameters
        
       userData.crtProc =  userData.procConstr(userData_main.MD, userData.crtPackage.outputDirectory_, ...
                                    userData_main.crtProc.maskRefineProcess_.funParams_);
    else
        % If segmentation process does not have a refinement process
        % Create a new refinement process with default parameters
        
        userData.crtProc = userData.procConstr(userData_main.MD, userData.crtPackage.outputDirectory_);
    end
    
    userData_main.maskRefineProc = userData.crtProc;
    
else
    userData.crtProc = userData_main.maskRefineProc;
end



% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;


% ---------------------- Parameter Setup -------------------------

funParams = userData.crtProc.funParams_;

if funParams.MaskCleanUp
    if ~funParams.FillHoles
        set(handles.checkbox_fillholes, 'Value', 0)
    end
    set(handles.edit_1, 'String',num2str(funParams.MinimumSize))
    set(handles.edit_2, 'String',num2str(funParams.ClosureRadius))
    set(handles.edit_3, 'String',num2str(funParams.ObjectNumber))
else
    set(handles.checkbox_cleanup, 'Value', 0)
    set(handles.checkbox_fillholes, 'Value', 0, 'Enable','off')
    set(handles.text_para1, 'Enable', 'off');
    set(handles.text_para2, 'Enable', 'off');
    set(handles.text_para3, 'Enable', 'off');
    set(handles.edit_1, 'Enable', 'off');
    set(handles.edit_2, 'Enable', 'off');
    set(handles.edit_3, 'Enable', 'off');
end

if funParams.EdgeRefinement
    set(handles.checkbox_edge, 'Value', 1)
    set(handles.text_para4, 'Enable', 'on');
    set(handles.text_para5, 'Enable', 'on');
    set(handles.text_para6, 'Enable', 'on');
    set(handles.edit_4, 'Enable', 'on', 'String',num2str(funParams.MaxEdgeAdjust));
    set(handles.edit_5, 'Enable', 'on', 'String',num2str(funParams.MaxEdgeGap));
    set(handles.edit_6, 'Enable', 'on', 'String',num2str(funParams.PreEdgeGrow));    
    
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
set(Img, 'UserData', userData.crtProc.getHelp(true))

% ----------------------------------------------------------------

% Update user data and GUI data
set(userData.mainFig, 'UserData', userData_main);
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = maskRefineProcessGUI_OutputFcn(hObject, eventdata, handles) 
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
userData = get(handles.figure1);
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button

userData = get(handles.figure1, 'UserData');

% -------- Check user input --------


if get(handles.checkbox_cleanup, 'value')
    if isnan(str2double(get(handles.edit_1, 'String'))) ...
            || str2double(get(handles.edit_1, 'String')) < 0
        errordlg('Please provide a valid input for ''Minimus Size''.','Setting Error','modal');
        return;
    end 
    if isnan(str2double(get(handles.edit_2, 'String'))) ...
            || str2double(get(handles.edit_2, 'String')) < 0
        errordlg('Please provide a valid input for ''Closure Radius''.','Setting Error','modal');
        return;
    end
    if isnan(str2double(get(handles.edit_3, 'String'))) ...
            || str2double(get(handles.edit_3, 'String')) < 0
        errordlg('Please provide a valid input for ''Object Number''.','Setting Error','modal');
        return;
    end     
end

if get(handles.checkbox_edge, 'value')
    if isnan(str2double(get(handles.edit_4, 'String'))) ...
            || str2double(get(handles.edit_4, 'String')) < 0
        errordlg('Please provide a valid input for ''Maximum Adjust Distance''.','Setting Error','modal');
        return;
    end 
    if isnan(str2double(get(handles.edit_5, 'String'))) ...
            || str2double(get(handles.edit_5, 'String')) < 0
        errordlg('Please provide a valid input for ''Maximum Edge Gap''.','Setting Error','modal');
        return;
    end
    if isnan(str2double(get(handles.edit_6, 'String'))) ...
            || str2double(get(handles.edit_6, 'String')) < 0
        errordlg('Please provide a valid input for ''Radius of Growth''.','Setting Error','modal');
        return;
    end     
end

if ~get(handles.checkbox_cleanup, 'value') && ~get(handles.checkbox_edge, 'value')
    errordlg('Please select at least one option for mask refinement processing.')
    return;
end
   
% -------- Set parameter --------

funParams = userData.crtProc.funParams_;
    
if userData.crtProc.procChanged_ 
    
    % Get parameter
    
    if get(handles.checkbox_cleanup, 'Value')
        funParams.MaskCleanUp = true;
        funParams.MinimumSize = str2double(get(handles.edit_1, 'String'));
        funParams.ClosureRadius = str2double(get(handles.edit_2, 'String'));
        funParams.ObjectNumber = str2double(get(handles.edit_3, 'String'));
        if get(handles.checkbox_fillholes, 'Value')
            funParams.FillHoles = true;
        else
            funParams.FillHoles = false;
        end
    else
        funParams.MaskCleanUp = false;
    end
    
    if get(handles.checkbox_edge, 'Value')
        funParams.EdgeRefinement = true;
        funParams.MaxEdgeAdjust = str2double(get(handles.edit_4, 'String'));
        funParams.MaxEdgeGap = str2double(get(handles.edit_5, 'String'));
        funParams.PreEdgeGrow = str2double(get(handles.edit_6, 'String'));
    else
        funParams.EdgeRefinement = false;
    end
       
    % Set parameters
    userData.crtProc.setPara(funParams);
end

% Save user data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);


% --- Executes on button press in checkbox_cleanup.
function checkbox_cleanup_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);

switch get(hObject, 'Value')
    case 0
        set(handles.text_para1, 'Enable', 'off');
        set(handles.edit_1,'Enable','off');
        set(handles.text_para2, 'Enable', 'off');
        set(handles.edit_2,'Enable','off');
        set(handles.text_para3, 'Enable', 'off');
        set(handles.edit_3,'Enable','off');
        set(handles.checkbox_fillholes,'Enable','off');          
    case 1
        set(handles.text_para1, 'Enable', 'on');
        set(handles.edit_1,'Enable','on');
        set(handles.text_para2, 'Enable', 'on');
        set(handles.edit_2,'Enable','on');
        set(handles.text_para3, 'Enable', 'on');
        set(handles.edit_3,'Enable','on');
        set(handles.checkbox_fillholes,'Enable','on');        
end



function edit_1_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes during object creation, after setting all properties.
function edit_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_2_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


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



function edit_3_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


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


% --- Executes on button press in checkbox_fillholes.
function checkbox_fillholes_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes on button press in checkbox_edge.
function checkbox_edge_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_edge (see GCBO)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);

switch get(hObject, 'Value')
    case 0
        set(handles.text_para4, 'Enable', 'off');
        set(handles.edit_4,'Enable','off');
        set(handles.text_para5, 'Enable', 'off');
        set(handles.edit_5,'Enable','off');
        set(handles.text_para6, 'Enable', 'off');
        set(handles.edit_6,'Enable','off');
        
    case 1
        set(handles.text_para4, 'Enable', 'on');
        set(handles.edit_4,'Enable','on');
        set(handles.text_para5, 'Enable', 'on');
        set(handles.edit_5,'Enable','on');
        set(handles.text_para6, 'Enable', 'on');
        set(handles.edit_6,'Enable','on');    
end



function edit_4_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


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



function edit_5_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes during object creation, after setting all properties.
function edit_5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_6_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData.crtProc.setProcChanged(true);


% --- Executes during object creation, after setting all properties.
function edit_6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
