function varargout = costMatRandomDirectedSwitchingMotionCloseGapsGUI(varargin)
% COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI M-file for costMatRandomDirectedSwitchingMotionCloseGapsGUI.fig
%      COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI, by itself, creates a new COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI or raises the existing
%      singleton*.
%
%      H = COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI returns the handle to a new COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI or the handle to
%      the existing singleton*.
%
%      COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI.M with the given input arguments.
%
%      COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI('Property','Value',...) creates a new COSTMATRANDOMDIRECTEDSWITCHINGMOTIONCLOSEGAPSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help costMatRandomDirectedSwitchingMotionCloseGapsGUI

% Last Modified by GUIDE v2.5 09-Dec-2011 16:37:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @costMatRandomDirectedSwitchingMotionCloseGapsGUI_OutputFcn, ...
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


% --- Executes just before costMatRandomDirectedSwitchingMotionCloseGapsGUI is made visible.
function costMatRandomDirectedSwitchingMotionCloseGapsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.gapclosingFig = costMatRandomDirectedSwitchingMotionCloseGapsGUI{procID}('mainFig',
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
set(handles.checkbox_linearMotion, 'Value', logical(parameters.linearMotion));
set(handles.checkbox_immediateDirectionReversal, 'Value', parameters.linearMotion==2);
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
    set(Img, 'UserData', struct('class', 'costMatLinearMotionCloseGaps2GUI'))
else
    set(Img, 'UserData', 'Please refer to help file.')
end



set(handles.figure1, 'UserData', userData)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes costMatRandomDirectedSwitchingMotionCloseGapsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = costMatRandomDirectedSwitchingMotionCloseGapsGUI_OutputFcn(hObject, eventdata, handles) 
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
