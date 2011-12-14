function varargout = costMatLinearMotionCloseGaps2GUI(varargin)
% COSTMATLINEARMOTIONCLOSEGAPS2GUI M-file for costMatLinearMotionCloseGaps2GUI.fig
%      COSTMATLINEARMOTIONCLOSEGAPS2GUI, by itself, creates a new COSTMATLINEARMOTIONCLOSEGAPS2GUI or raises the existing
%      singleton*.
%
%      H = COSTMATLINEARMOTIONCLOSEGAPS2GUI returns the handle to a new COSTMATLINEARMOTIONCLOSEGAPS2GUI or the handle to
%      the existing singleton*.
%
%      COSTMATLINEARMOTIONCLOSEGAPS2GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COSTMATLINEARMOTIONCLOSEGAPS2GUI.M with the given input arguments.
%
%      COSTMATLINEARMOTIONCLOSEGAPS2GUI('Property','Value',...) creates a new COSTMATLINEARMOTIONCLOSEGAPS2GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before costMatLinearMotionCloseGaps2GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to costMatLinearMotionCloseGaps2GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help costMatLinearMotionCloseGaps2GUI

% Last Modified by GUIDE v2.5 21-Mar-2011 17:13:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @costMatLinearMotionCloseGaps2GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @costMatLinearMotionCloseGaps2GUI_OutputFcn, ...
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


% --- Executes just before costMatLinearMotionCloseGaps2GUI is made visible.
function costMatLinearMotionCloseGaps2GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.gapclosingFig = costMatLinearMotionCloseGaps2GUI{procID}('mainFig',
% handles.figure1, procID);
%
% userData.mainFig
% userData.procID
% userData.handles_main
% userData.userData_main
% userData.crtProc
% userData.parameters

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

handles.output = hObject;
userData = get(handles.figure1, 'UserData');

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);
userData.userData_main = get(userData.handles_main.figure1, 'UserData');
userData.crtProc = userData.userData_main.crtProc;

u = get(userData.handles_main.popupmenu_gapclosing, 'UserData');
userData.parameters = u{userData.procID};
parameters = userData.parameters;

% Parameter Setup
set(handles.checkbox_linearMotion, 'Value', parameters.linearMotion)
if parameters.linearMotion
    
    arrayfun(@(x)eval(['set(handles.text_linear_',num2str(x),', ''Enable'', ''on'')']), 1:7)
    set(handles.edit_lenForClassify , 'Enable', 'on')
    set(handles.edit_linStdMult , 'Enable', 'on')
    set(handles.edit_gapLengthTransitionL , 'Enable', 'on')
    set(handles.edit_before_2 , 'Enable', 'on')
    set(handles.edit_after_2 , 'Enable', 'on')
    set(handles.edit_maxAngleVV , 'Enable', 'on')
    
    set(handles.edit_lenForClassify, 'String', num2str(parameters.lenForClassify))
    set(handles.edit_linStdMult, 'String', num2str(parameters.linStdMult(1)))
    set(handles.edit_before_2, 'String', num2str(parameters.linScaling(1)))
    set(handles.edit_after_2, 'String', num2str(parameters.linScaling(2)))
    set(handles.edit_gapLengthTransitionL, 'String', num2str(parameters.timeReachConfL-1))  
    set(handles.edit_maxAngleVV, 'String', num2str(parameters.maxAngleVV))    
end

set(handles.edit_lower, 'String', num2str(parameters.minSearchRadius))
set(handles.edit_upper, 'String', num2str(parameters.maxSearchRadius))
set(handles.edit_brownStdMult, 'String', num2str(parameters.brownStdMult(1)))
set(handles.checkbox_useLocalDensity, 'Value', parameters.useLocalDensity)
set(handles.edit_nnWindow, 'String', num2str(parameters.nnWindow))
set(handles.edit_before, 'String', num2str(parameters.brownScaling(1)))
set(handles.edit_after, 'String', num2str(parameters.brownScaling(2)))
set(handles.edit_gapLengthTransitionB, 'String', num2str(parameters.timeReachConfB-1))

if isempty(parameters.ampRatioLimit) || (length(parameters.ampRatioLimit) ==1 && parameters.ampRatioLimit == 0)
    
    set(handles.checkbox_ampRatioLimit, 'Value', 0)
    arrayfun(@(x)eval( ['set(handles.text_ampRatioLimit_',num2str(x),', ''Enable'', ''off'')'] ), 1:3)
    set(handles.edit_min, 'Enable', 'off')
    set(handles.edit_max, 'Enable', 'off')
else
    
    set(handles.edit_min, 'String', num2str(parameters.ampRatioLimit(1)))
    set(handles.edit_max, 'String', num2str(parameters.ampRatioLimit(2)))
end

set(handles.edit_resLimit, 'String', num2str(parameters.resLimit))
set(handles.edit_gapPenalty, 'String', num2str(parameters.gapPenalty))

% Get icon infomation
userData.questIconData = userData.userData_main.questIconData;
userData.colormap = userData.userData_main.colormap;

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
    set(Img, 'UserData', struct('class', mfilename))
else
    set(Img, 'UserData', 'Please refer to help file.')
end



set(handles.figure1, 'UserData', userData)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes costMatLinearMotionCloseGaps2GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = costMatLinearMotionCloseGaps2GUI_OutputFcn(hObject, eventdata, handles) 
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
delete(handles.figure1)

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');
parameters = userData.parameters;

% if userData.crtProc.procChanged_ 
    
    minSearchRadius = get(handles.edit_lower, 'String');
    maxSearchRadius = get(handles.edit_upper, 'String');
    brownStdMult = get(handles.edit_brownStdMult, 'String');
    nnWindow = get(handles.edit_nnWindow, 'String');
    brownScaling_1 = get(handles.edit_before, 'String'); 
    brownScaling_2 = get(handles.edit_after, 'String'); 
    gapLengthTransitionB = get(handles.edit_gapLengthTransitionB, 'String'); 
    ampRatioLimit_1 = get(handles.edit_min, 'String'); 
    ampRatioLimit_2 = get(handles.edit_max, 'String'); 
    resLimit = get(handles.edit_resLimit, 'String'); 
    gapPenalty = get(handles.edit_gapPenalty, 'String'); 
    
    lenForClassify = get(handles.edit_lenForClassify, 'String'); 
    linStdMult = get(handles.edit_linStdMult, 'String'); 
    linScaling_1 = get(handles.edit_before_2, 'String'); 
    linScaling_2 = get(handles.edit_after_2, 'String');    
    gapLengthTransitionL = get(handles.edit_gapLengthTransitionL, 'String'); 
    maxAngleVV = get(handles.edit_maxAngleVV, 'String'); 
    
    % lower
    if isempty( minSearchRadius )
        errordlg('Parameter "Lower Bound" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(minSearchRadius)) || str2double(minSearchRadius) < 0
        errordlg('Please provide a valid value to parameter "Lower Bound".','Error','modal')
        return
    else
        minSearchRadius = str2double(minSearchRadius);
    end      
    
    % Upper
    if isempty( maxSearchRadius )
        errordlg('Parameter "Upper Bound" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(maxSearchRadius)) || str2double(maxSearchRadius) < 0 
        errordlg('Please provide a valid value to parameter "Upper Bound".','Error','modal')
        return
        
    elseif str2double(maxSearchRadius) < minSearchRadius
        errordlg('"Upper Bound" should be larger than "Lower Bound".','Error','modal')
        return
        
    else
        maxSearchRadius = str2double(maxSearchRadius);
    end        
    
    % brownStdMult
    if isempty( brownStdMult )
        errordlg('Parameter "Multiplication Factor for Search Radius Calculation" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(brownStdMult)) || str2double(brownStdMult) < 0
        errordlg('Please provide a valid value to parameter "Multiplication Factor for Search Radius Calculation".','Error','modal')
        return
    else
        brownStdMult = str2double(brownStdMult)*ones(userData.crtProc.funParams_.gapCloseParam.timeWindow,1);
    end  
    
    % nnWindow
    if isempty( nnWindow )
        errordlg('Parameter "Number of Frames for Nearest Neighbor Distance Calculation" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(nnWindow)) || str2double(nnWindow) < 0
        errordlg('Please provide a valid value to parameter "Number of Frames for Nearest Neighbor Distance Calculation".','Error','modal')
        return
    else
        nnWindow = str2double(nnWindow);
    end    
    
    % brownScaling
    if isempty( brownScaling_1 )
        errordlg('Parameter "Scaling Power in Fast Expansion Phase" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(brownScaling_1)) || str2double(brownScaling_1) < 0
        errordlg('Please provide a valid value to parameter "Scaling Power in Fast Expansion Phase".','Error','modal')
        return
    else
        brownScaling_1 = str2double(brownScaling_1);
    end        
    
    % brownScaling
    if isempty( brownScaling_2 )
        errordlg('Parameter "Scaling Power in Slow Expansion Phase" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(brownScaling_2)) || str2double(brownScaling_2) < 0
        errordlg('Please provide a valid value to parameter "Scaling Power in Slow Expansion Phase".','Error','modal')
        return
    else
        brownScaling_2 = str2double(brownScaling_2);
    end   
    
    brownScaling = [brownScaling_1 brownScaling_2];
    
    % gapLengthTransitionB
    if isempty( gapLengthTransitionB )
        errordlg('Parameter "Gap length to transition from Fast to Slow Expansion" is requied by the algorithm.','Error','modal')
        return

    elseif isnan(str2double(gapLengthTransitionB)) || str2double(gapLengthTransitionB) < 0
        errordlg('Please provide a valid value to parameter "Gap length to transition from Fast to Slow Expansion".','Error','modal')
        return
    else
        gapLengthTransitionB = str2double(gapLengthTransitionB);
    end      
     
    % ampRatioLimit
    if ~get(handles.checkbox_ampRatioLimit, 'Value')
        ampRatioLimit = [];
    else
        % ampRatioLimit_1
        if isempty( ampRatioLimit_1 )
            errordlg('Parameter "Min Allowed" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(ampRatioLimit_1)) || str2double(ampRatioLimit_1) < 0
            errordlg('Please provide a valid value to parameter "Min Allowed".','Error','modal')
            return
        else
            ampRatioLimit_1 = str2double(ampRatioLimit_1);
        end        

        % ampRatioLimit_2
        if isempty( ampRatioLimit_2 )
            errordlg('Parameter "Max Allowed" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(ampRatioLimit_2)) || str2double(ampRatioLimit_2) < 0
            errordlg('Please provide a valid value to parameter "Max Allowed".','Error','modal')
            return
            
        elseif str2double(ampRatioLimit_2) <= ampRatioLimit_1
            errordlg('"Max Allowed" should be larger than "Min Allowed".','Error','modal')
            return   
            
        else
            ampRatioLimit_2 = str2double(ampRatioLimit_2);
        end   

        ampRatioLimit = [ampRatioLimit_1 ampRatioLimit_2];        

    end
    
    % resLimit
    if isempty( resLimit )
        resLimit = [];

    elseif isnan(str2double(resLimit)) || str2double(resLimit) < 0
        errordlg('Please provide a valid value to parameter "Time to Reach Confinement".','Error','modal')
        return
    else
        resLimit = str2double(resLimit);
    end      
    
    % gapPenalty
    if isempty( gapPenalty )
        gapPenalty = [];

    elseif isnan(str2double(gapPenalty)) || str2double(gapPenalty) < 0
        errordlg('Please provide a valid value to parameter "Time to Reach Confinement".','Error','modal')
        return
    else
        gapPenalty = str2double(gapPenalty);
    end    
    
    % If parameters.linearMotion = 1
    if parameters.linearMotion
       
        % lenForClassify
        if isempty( lenForClassify )
            errordlg('Parameter "Minimum Track Segment Length to Classify it as Linear or Random" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(lenForClassify)) || str2double(lenForClassify) < 0
            errordlg('Please provide a valid value to parameter "Minimum Track Segment Length to Classify it as Linear or Random".','Error','modal')
            return
        else
            lenForClassify = str2double(lenForClassify);
        end    
        
        % linStdMult
        if isempty( linStdMult )
            errordlg('Parameter "Multiplication Factor for Linear Search Radius Calculation" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(linStdMult)) || str2double(linStdMult) < 0
            errordlg('Please provide a valid value to parameter "Multiplication Factor for Linear Search Radius Calculation".','Error','modal')
            return
        else
            linStdMult = str2double(linStdMult)*ones(userData.crtProc.funParams_.gapCloseParam.timeWindow,1);
        end          
        
        % linScaling_1
        if isempty( linScaling_1 )
            errordlg('Parameter "Scaling Power in Fast Expansion Phase" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(linScaling_1)) || str2double(linScaling_1) < 0
            errordlg('Please provide a valid value to parameter "Scaling Power in Fast Expansion Phase".','Error','modal')
            return
        else
            linScaling_1 = str2double(linScaling_1);
        end        

        % linScaling_1
        if isempty( linScaling_2 )
            errordlg('Parameter "Scaling Power in Slow Expansion Phase" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(linScaling_2)) || str2double(linScaling_2) < 0
            errordlg('Please provide a valid value to parameter "Scaling Power in Slow Expansion Phase".','Error','modal')
            return
        else
            linScaling_2 = str2double(linScaling_2);
        end   

        linScaling = [linScaling_1 linScaling_2];  
        
        % gapLengthTransitionL
        if isempty( gapLengthTransitionL )
            errordlg('Parameter "Gap length to transition from Fast to Slow Expansion" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(gapLengthTransitionL)) || str2double(gapLengthTransitionL) < 0
            errordlg('Please provide a valid value to parameter "Gap length to transition from Fast to Slow Expansion".','Error','modal')
            return
        else
            gapLengthTransitionL = str2double(gapLengthTransitionL);
        end      
        
        % maxAngleVV
        if isempty( maxAngleVV )
            errordlg('Parameter "Maximum Angle Between Linear Track Segments" is requied by the algorithm.','Error','modal')
            return

        elseif isnan(str2double(maxAngleVV)) || str2double(maxAngleVV) < 0
            errordlg('Please provide a valid value to parameter "Maximum Angle Between Linear Track Segments".','Error','modal')
            return
        else
            maxAngleVV = str2double(maxAngleVV);
        end           
        
    end
    
    % ----------- Set Parameters --------------
    
    parameters.minSearchRadius = minSearchRadius;
    parameters.maxSearchRadius = maxSearchRadius;
    parameters.brownStdMult = brownStdMult;
    parameters.useLocalDensity = get(handles.checkbox_useLocalDensity, 'Value');
    parameters.nnWindow = nnWindow;
    parameters.brownScaling = brownScaling;
    parameters.timeReachConfB = gapLengthTransitionB+1;
    parameters.ampRatioLimit = ampRatioLimit;
    parameters.resLimit = resLimit;
    parameters.gapPenalty = gapPenalty;
    
    if parameters.linearMotion
        
        parameters.lenForClassify = lenForClassify;
        parameters.linStdMult = linStdMult;
        parameters.linScaling = linScaling;
        parameters.timeReachConfL = gapLengthTransitionL+1;
        parameters.maxAngleVV = maxAngleVV;
    end
    
    u = get(userData.handles_main.popupmenu_gapclosing, 'UserData');
    u{userData.procID} = parameters;
    
    set(userData.handles_main.popupmenu_gapclosing, 'UserData', u)    
    
% end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_brownStdMult_Callback(hObject, eventdata, handles)
% hObject    handle to edit_brownStdMult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_brownStdMult as text
%        str2double(get(hObject,'String')) returns contents of edit_brownStdMult as a double


% --- Executes during object creation, after setting all properties.
function edit_brownStdMult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_brownStdMult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function edit_lenForClassify_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lenForClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lenForClassify as text
%        str2double(get(hObject,'String')) returns contents of edit_lenForClassify as a double


% --- Executes during object creation, after setting all properties.
function edit_lenForClassify_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lenForClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_linStdMult_Callback(hObject, eventdata, handles)
% hObject    handle to edit_linStdMult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_linStdMult as text
%        str2double(get(hObject,'String')) returns contents of edit_linStdMult as a double


% --- Executes during object creation, after setting all properties.
function edit_linStdMult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_linStdMult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_maxAngleVV_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxAngleVV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxAngleVV as text
%        str2double(get(hObject,'String')) returns contents of edit_maxAngleVV as a double


% --- Executes during object creation, after setting all properties.
function edit_maxAngleVV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxAngleVV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_before_Callback(hObject, eventdata, handles)
% hObject    handle to edit_before (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_before as text
%        str2double(get(hObject,'String')) returns contents of edit_before as a double


% --- Executes during object creation, after setting all properties.
function edit_before_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_before (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_after_Callback(hObject, eventdata, handles)
% hObject    handle to edit_after (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_after as text
%        str2double(get(hObject,'String')) returns contents of edit_after as a double


% --- Executes during object creation, after setting all properties.
function edit_after_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_after (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gapLengthTransitionB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gapLengthTransitionB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gapLengthTransitionB as text
%        str2double(get(hObject,'String')) returns contents of edit_gapLengthTransitionB as a double


% --- Executes during object creation, after setting all properties.
function edit_gapLengthTransitionB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gapLengthTransitionB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ampRatioLimit.
function checkbox_ampRatioLimit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ampRatioLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ampRatioLimit

if get(hObject, 'Value')
    arrayfun(@(x)eval( ['set(handles.text_ampRatioLimit_',num2str(x),', ''Enable'', ''on'')'] ), 1:3)
    set(handles.edit_min, 'Enable', 'on')
    set(handles.edit_max, 'Enable', 'on')    
else
    arrayfun(@(x)eval( ['set(handles.text_ampRatioLimit_',num2str(x),', ''Enable'', ''off'')'] ), 1:3)
    set(handles.edit_min, 'Enable', 'off')
    set(handles.edit_max, 'Enable', 'off')    
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gapPenalty_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gapPenalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gapPenalty as text
%        str2double(get(hObject,'String')) returns contents of edit_gapPenalty as a double


% --- Executes during object creation, after setting all properties.
function edit_gapPenalty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gapPenalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_resLimit_Callback(hObject, eventdata, handles)
% hObject    handle to edit_resLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_resLimit as text
%        str2double(get(hObject,'String')) returns contents of edit_resLimit as a double


% --- Executes during object creation, after setting all properties.
function edit_resLimit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_resLimit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox_linearMotion.
function checkbox_linearMotion_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_linearMotion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_linearMotion



function edit_lower_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lower as text
%        str2double(get(hObject,'String')) returns contents of edit_lower as a double


% --- Executes during object creation, after setting all properties.
function edit_lower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_upper_Callback(hObject, eventdata, handles)
% hObject    handle to edit_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_upper as text
%        str2double(get(hObject,'String')) returns contents of edit_upper as a double


% --- Executes during object creation, after setting all properties.
function edit_upper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_useLocalDensity.
function checkbox_useLocalDensity_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_useLocalDensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_useLocalDensity



function edit_nnWindow_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nnWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nnWindow as text
%        str2double(get(hObject,'String')) returns contents of edit_nnWindow as a double


% --- Executes during object creation, after setting all properties.
function edit_nnWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nnWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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



function edit_before_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_before_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_before_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_before_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_before_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_before_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_after_2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_after_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_after_2 as text
%        str2double(get(hObject,'String')) returns contents of edit_after_2 as a double


% --- Executes during object creation, after setting all properties.
function edit_after_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_after_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gapLengthTransitionL_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gapLengthTransitionL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gapLengthTransitionL as text
%        str2double(get(hObject,'String')) returns contents of edit_gapLengthTransitionL as a double


% --- Executes during object creation, after setting all properties.
function edit_gapLengthTransitionL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gapLengthTransitionL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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
