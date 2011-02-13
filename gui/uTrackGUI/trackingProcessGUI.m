function varargout = trackingProcessGUI(varargin)
% TRACKINGPROCESSGUI M-file for trackingProcessGUI.fig
%      TRACKINGPROCESSGUI, by itself, creates a new TRACKINGPROCESSGUI or raises the existing
%      singleton*.
%
%      H = TRACKINGPROCESSGUI returns the handle to a new TRACKINGPROCESSGUI or the handle to
%      the existing singleton*.
%
%      TRACKINGPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKINGPROCESSGUI.M with the given input arguments.
%
%      TRACKINGPROCESSGUI('Property','Value',...) creates a new TRACKINGPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackingProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackingProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackingProcessGUI

% Last Modified by GUIDE v2.5 15-Dec-2010 15:20:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackingProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @trackingProcessGUI_OutputFcn, ...
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


% --- Executes just before trackingProcessGUI is made visible.
function trackingProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% userData.setFig(procID) = trackingProcessGUI('mainFig',handles.figure1, procID);
%
% Available tools 
% UserData data:
%       userData.MD - 1x1 the current movie data
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main figure
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%
%       userData.cost_linking
%       userData.cost_gapclosing
%       userData.kalman_reserveMem
%       userData.kalman_initialize
%       userData.kalman_calcGain
%       userData.kalman_timeReverse
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig - the handle of setting panel for kalman filter initilization
%
%       userData.path
%       userData.file

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
userData.MD = userData_main.MD(userData_main.id);  % Get the current Movie Data
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get current process constructer
eval ( [ 'userData.procConstr = @', ...
    userData.crtPackage.processClassNames_{userData.procID},';']);

% Set cost function and kalman options
userData.cost_linking = {'costMatLinearMotionLink2'};
userData.cost_gapclosing = {'costMatLinearMotionCloseGaps2'};

userData.fun_cost_linking = {@costMatLinearMotionLink2GUI};
userData.fun_cost_gap = {@costMatLinearMotionCloseGaps2GUI};

userData.kalman_reserveMem = {'kalmanResMemLM'};
userData.kalman_initialize = {'kalmanInitLinearMotion'};
userData.fun_kalman_initialize = {@kalmanInitializationGUI};
userData.kalman_calcGain = {'kalmanGainLinearMotion'};
userData.kalman_timeReverse = {'kalmanReverseLinearMotion'};

% If process does not exist, create a default one in user data.
if isempty(userData.crtProc)
    userData.crtProc = userData.procConstr(userData.MD, userData.crtPackage.outputDirectory_);                       
end

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% ---------------------- Channel Setup -------------------------

channelIndex = userData.crtProc.channelIndex_;
% Set up available input channels
set(handles.listbox_1, 'String', {userData.MD.channels_.channelPath_},...
        'Userdata', 1: length(userData.MD.channels_));
    
set(handles.listbox_2, 'String', ...
        {userData.MD.channels_(channelIndex).channelPath_}, ...
        'Userdata',channelIndex);    

% ---------------------- Parameter Setup -------------------------

funParams = userData.crtProc.funParams_;

set(handles.edit_probDim, 'String', num2str(funParams.probDim))
set(handles.checkbox_verbose, 'Value', funParams.verbose)

% gapCloseParam
set(handles.edit_maxgap, 'String', num2str(funParams.gapCloseParam.timeWindow - 1))
set(handles.edit_minlength, 'String', num2str(funParams.gapCloseParam.minTrackLen))
set(handles.checkbox_histogram, 'Value', funParams.gapCloseParam.diagnostics)
set(handles.checkbox_overwrite, 'Value', userData.crtProc.overwrite_)

if funParams.gapCloseParam.mergeSplit == 1
    set(handles.checkbox_merging, 'Value',1), set(handles.checkbox_splitting, 'Value',1)
    
elseif funParams.gapCloseParam.mergeSplit == 2
    set(handles.checkbox_merging, 'Value',1), set(handles.checkbox_splitting, 'Value',0)
    
elseif funParams.gapCloseParam.mergeSplit == 3
    set(handles.checkbox_merging, 'Value',0), set(handles.checkbox_splitting, 'Value',1)
    
elseif funParams.gapCloseParam.mergeSplit == 0
    set(handles.checkbox_merging, 'Value',0), set(handles.checkbox_splitting, 'Value',0)
    
else
   error('User-defined: Wrong parameter.') 
end
    
% costMatrices
i1 = find(strcmp(funParams.costMatrices(1).funcName, userData.cost_linking));
i2 = find(strcmp(funParams.costMatrices(2).funcName, userData.cost_gapclosing));

if length(i1)>1 || length(i2)>1 || isempty(i1) || isempty(i2)
    error('User-defined: the length of matching methods must be 1.')
end

u1 = cell(1, length(get(handles.popupmenu_linking, 'String')));
u2 = cell(1, length(get(handles.popupmenu_gapclosing, 'String')));

u1{i1} = funParams.costMatrices(1).parameters;
u2{i2} = funParams.costMatrices(2).parameters;


set(handles.popupmenu_linking, 'Value', i1, 'UserData', u1)
set(handles.popupmenu_gapclosing, 'Value', i2, 'UserData', u2)


% kalmanFunctions
i1 = find(strcmp(funParams.kalmanFunctions.reserveMem, userData.kalman_reserveMem));
i2 = find(strcmp(funParams.kalmanFunctions.initialize, userData.kalman_initialize));
i3 = find(strcmp(funParams.kalmanFunctions.calcGain, userData.kalman_calcGain));
i4 = find(strcmp(funParams.kalmanFunctions.timeReverse, userData.kalman_timeReverse));

if length(i1)>1 || length(i2)>1 || isempty(i1) || isempty(i2) || length(i3)>1 || length(i4)>1 || isempty(i3) || isempty(i4)
    error('User-defined: the length of matching methods must be 1.')
end

u2 = cell(1, length(get(handles.popupmenu_kalman_initialize, 'String')));
u2{i2} = funParams.costMatrices(1).parameters.kalmanInitParam;

set(handles.popupmenu_kalman_reserve, 'Value', i1)
set(handles.popupmenu_kalman_initialize, 'Value', i2, 'UserData', u2)
set(handles.popupmenu_kalman_gain, 'Value', i3)
set(handles.popupmenu_kalman_reverse, 'Value', i4)

% funParams.saveResults

userData.path = funParams.saveResults.dir;
userData.file = userData.crtProc.filename_;

% Show the actual file name on GUI
i = userData.crtProc.channelIndex_(1);
file = ['Channel_' num2str(i) '_' userData.crtProc.filename_];
if ~userData.crtProc.overwrite_
    file = enumFileName(userData.path, file);
end

str = [funParams.saveResults.dir file];
if length(str)>100
    str = ['...' str(end-100:end)];
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
    set(Img, 'UserData', struct('class', 'TrackingProcess'))
else
    set(Img, 'UserData', userData.crtProc.getHelp(true))
end


% ----------------------------------------------------------------

% Update user data and GUI data
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = trackingProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
userData_main = get(userData.mainFig, 'UserData');

% --------------- Check User Input ---------------

if isempty(get(handles.listbox_2, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
end

probDim = get(handles.edit_probDim, 'String');
timeWindow = get(handles.edit_maxgap, 'String');
minTrackLen = get(handles.edit_minlength, 'String');

if isempty( probDim )
    errordlg('Parameter "Problem Dimensionality" is required by the algorithm.','Error','modal')
    return
    
elseif isnan(str2double(probDim)) || str2double(probDim) < 0 || floor(str2double(probDim)) ~= ceil(str2double(probDim))
    errordlg('Please provide a valid value to parameter "Problem Dimensionality".','Error','modal')
    return
    
else
    probDim = str2double(probDim);
end

if isempty( timeWindow )
    errordlg('Parameter "Maximum Gap to Close" is required by the algorithm.','Error','modal')
    return
    
elseif isnan(str2double(timeWindow)) || str2double(timeWindow) < 0 || floor(str2double(timeWindow)) ~= ceil(str2double(timeWindow))
    errordlg('Please provide a valid value to parameter "Maximum Gap to Close".','Error','modal')
    return
    
else
    timeWindow = str2double(timeWindow) + 1;
end

if isempty( minTrackLen )
    errordlg('Parameter "Minimum Length of Track Segment from First Step to use in Second Step" is required by the algorithm.','Error','modal')
    return
    
elseif isnan(str2double(minTrackLen)) || str2double(minTrackLen) < 0 || floor(str2double(minTrackLen)) ~= ceil(str2double(minTrackLen))
    errordlg('Please provide a valid value to parameter "Minimum Length of Track Segment from First Step to use in Second Step".','Error','modal')
    return
    
else
    minTrackLen = str2double(minTrackLen);
end

i_linking = get(handles.popupmenu_linking, 'Value');
i_gapclosing = get(handles.popupmenu_gapclosing, 'Value');
i_kalman = get(handles.popupmenu_kalman_initialize, 'Value');

u_linking = get(handles.popupmenu_linking, 'UserData');
u_gapclosing = get(handles.popupmenu_gapclosing, 'UserData');
u_kalman = get(handles.popupmenu_kalman_initialize, 'UserData');

if isempty( u_linking{i_linking} )
    
    errordlg('Plese set up the selected cost function for "Step 1: frame-to-frame linking".','Error','modal')
end

if isempty( u_gapclosing{i_gapclosing} )
    
    errordlg('Plese set up the selected cost function for "Step 2: gap closing, mergin and splitting".','Error','modal')
end

%---------Check if channel indexs are changed---------

channelIndex = get (handles.listbox_2, 'Userdata');

if ~isempty( setdiff(channelIndex, userData.crtProc.channelIndex_) ) ...
    || ~isempty( setdiff(userData.crtProc.channelIndex_, channelIndex) )

    % If channel indexs are changed, set procChanged to true
    userData.crtProc.setProcChanged(true);
end

% -------- Set parameter --------

funParams = userData.crtProc.funParams_;

% Set channels
userData.crtProc.setChannelIndex(channelIndex)

funParams.probDim = probDim;
funParams.verbose = get(handles.checkbox_verbose, 'Value');
funParams.gapCloseParam.timeWindow = timeWindow;
funParams.gapCloseParam.minTrackLen = minTrackLen;
funParams.gapCloseParam.diagnostics = get(handles.checkbox_histogram, 'Value');

if get(handles.checkbox_merging, 'Value') && get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 1;
elseif get(handles.checkbox_merging, 'Value') && ~get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 2;
elseif ~get(handles.checkbox_merging, 'Value') && get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 3;
elseif ~get(handles.checkbox_merging, 'Value') && ~get(handles.checkbox_splitting, 'Value')
    funParams.gapCloseParam.mergeSplit = 0;
end

funParams.saveResults.dir = userData.path;

funParams.costMatrices(1).funcName = userData.cost_linking{i_linking};
funParams.costMatrices(1).parameters = u_linking{i_linking};
funParams.costMatrices(2).funcName = userData.cost_gapclosing{i_gapclosing};
funParams.costMatrices(2).parameters = u_gapclosing{i_gapclosing};

funParams.kalmanFunctions.initialize = userData.kalman_initialize{i_kalman};
funParams.costMatrices(1).parameters.kalmanInitParam = u_kalman{i_kalman};

% Set up parameters effected by funParams.gapCloseParam.timeWindow
funParams.costMatrices(2).parameters.brownStdMult = funParams.costMatrices(2).parameters.brownStdMult(1) * ones(funParams.gapCloseParam.timeWindow,1);
funParams.costMatrices(2).parameters.linStdMult = funParams.costMatrices(2).parameters.linStdMult(1) * ones(funParams.gapCloseParam.timeWindow,1);

% Save result file name
userData.crtProc.setFileName(userData.file)

% Set Overwrite flag
userData.crtProc.setOverwrite(get(handles.checkbox_overwrite, 'Value'))

% Save funParams
userData.crtProc.setPara(funParams);

% -------------------------- Assign process ----------------------------

% If this is a brand new process, attach current process to MovieData and 
% package's process list 
if isempty( userData.crtPackage.processes_{userData.procID} )
    
    % Add new process to both process lists of MovieData and current package
    userData_main.MD(userData_main.id).addProcess( userData.crtProc );
    userData.crtPackage.setProcess(userData.procID, userData.crtProc);
    
    % Set font weight of process name bold
    eval([ 'set(userData.handles_main.checkbox_',...
            num2str(userData.procID),', ''FontWeight'',''bold'')' ]);
end

% ----------------------Sanity Check (II, III check)----------------------

% Do sanity check - only check changed parameters
procEx = userData.crtPackage.sanityCheck(false,'all');

% Return user data !!!
set(userData.mainFig, 'UserData', userData_main)

% Draw some bugs on the wall 
for i = 1: length(procEx)
   if ~isempty(procEx{i})
       % Draw warning label on the i th process
       userfcn_drawIcon(userData.handles_main,'warn',i,procEx{i}(1).message, true) % user data is retrieved, updated and submitted
   end
end
% Refresh user data !!
userData_main = get(userData.mainFig, 'UserData');


% -------------------- Apply setting to all movies ------------------------

if get(handles.checkbox_applytoall, 'Value')
    
for x = 1: length(userData_main.MD)
   
   if x == userData_main.id
      continue 
   end
   
   % set segmentation process' parameters:
   % ChannelIndex - all channels
   % funParams.saveResults.dir - result directory
   % userData.crtProc.filename_ result file name
   % overwrite_
   
   % Channel Index
   l = length(userData_main.MD(x).channels_);
   temp = arrayfun(@(x)(x > l),channelIndex, 'UniformOutput', true );
   channelIndex_ = channelIndex(logical(~temp));   

   % output dir
   funParams.saveResults.dir = [userData_main.package(x).outputDirectory_  filesep userData.crtProc.name_ filesep];
   
   % if new process, create a new process with funParas and add to
   % MovieData and package's process list
   if isempty(userData_main.package(x).processes_{userData.procID})
       
       process = userData.procConstr(userData_main.MD(x), userData_main.package(x).outputDirectory_, channelIndex_, funParams);
       userData_main.MD(x).addProcess( process )
       userData_main.package(x).setProcess(userData.procID, process )
   else
       userData_main.package(x).processes_{userData.procID}.setPara(funParams)
       userData_main.package(x).processes_{userData.procID}.setChannelIndex(channelIndex_)
   end
   
   userData_main.package(x).processes_{userData.procID}.setFileName(userData.crtProc.filename_)
   userData_main.package(x).processes_{userData.procID}.setOverwrite(get(handles.checkbox_overwrite, 'Value'))
   
   % If current process is changed, then assume funParams are changed in
   % all movies
   if userData.crtProc.procChanged_ 
       
       userData_main.package(x).processes_{userData.procID}.setProcChanged(true);
   end
   
    % Do sanity check - only check changed parameters
    procEx = userData_main.package(x).sanityCheck(false,'all');

    % Record the exceptions
    for i = 1: length(procEx)
       if ~isempty(procEx{i})
           % Record the icon and message to user data
           userData_main.statusM(x).IconType{userData.procID} = 'warn';
           userData_main.statusM(x).Msg{i} = procEx{i}(1).message;
       end
    end      
    
end
    
% Save user data
set(userData.mainFig, 'UserData', userData_main)    

end

% Save user data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userData = get(handles.figure1, 'UserData');

fig(1) = isfield(userData, 'linkingFig') && ishandle(userData.linkingFig);
fig(2) = isfield(userData, 'gapclosingFig') && ishandle(userData.gapclosingFig);
fig(3) = isfield(userData, 'kalmanFig') && ishandle(userData.kalmanFig);

if any(fig)
    
    temp = find(fig);
    
    switch temp(1)
        
        case 1
            hPopupmenu = handles.popupmenu_linking;
        case 2
            hPopupmenu = handles.popupmenu_gapclosing;
        case 3
            hPopupmenu = handles.popupmenu_kalman_initialize;
    end
    str = get(hPopupmenu, 'String');
    user_response = questdlg(['The setting panel of method "',str{get(hPopupmenu, 'Value')},'" is still open. Do you want to close the setting panel of tracking?'], ...
        'Setting Open', 'No','Yes','Yes');
    
    if strcmpi(user_response, 'yes')

        delete(handles.figure1)
    end   
    
else

    delete(handles.figure1);
end

% --- Executes on button press in checkbox_applytoall.
function checkbox_applytoall_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_applytoall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_applytoall


% --- Executes on selection change in popupmenu_kalman_reserve.
function popupmenu_kalman_reserve_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_reserve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_kalman_reserve contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_kalman_reserve


% --- Executes during object creation, after setting all properties.
function popupmenu_kalman_reserve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_reserve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_kalman_initialize.
function popupmenu_kalman_initialize_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_kalman_initialize contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_kalman_initialize


% --- Executes during object creation, after setting all properties.
function popupmenu_kalman_initialize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_initialize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_kalman_gain.
function popupmenu_kalman_gain_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_kalman_gain contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_kalman_gain


% --- Executes during object creation, after setting all properties.
function popupmenu_kalman_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_kalman_reverse.
function popupmenu_kalman_reverse_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_reverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_kalman_reverse contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_kalman_reverse


% --- Executes during object creation, after setting all properties.
function popupmenu_kalman_reverse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_kalman_reverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_linking.
function popupmenu_linking_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_linking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_linking contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_linking


% --- Executes during object creation, after setting all properties.
function popupmenu_linking_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_linking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_set_linking.
function pushbutton_set_linking_Callback(hObject, eventdata, handles)
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig

userData = get(handles.figure1, 'UserData');
procID = get(handles.popupmenu_linking, 'Value');
if procID > length(userData.fun_cost_linking)
    warndlg('Please select a cost function for linking step.','Error','modal')
    return
else
    userData.linkingFig = userData.fun_cost_linking{procID}('mainFig', handles.figure1, procID);
end
set(handles.figure1, 'UserData', userData);

% --- Executes on selection change in popupmenu_gapclosing.
function popupmenu_gapclosing_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_gapclosing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_gapclosing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_gapclosing


% --- Executes during object creation, after setting all properties.
function popupmenu_gapclosing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_gapclosing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_set_gapclosing.
function pushbutton_set_gapclosing_Callback(hObject, eventdata, handles)
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig
userData = get(handles.figure1, 'UserData');
procID = get(handles.popupmenu_gapclosing, 'Value');
if procID > length(userData.fun_cost_gap)
    warndlg('Please select a cost function for gap closing step.','Error','modal')
    return
else
    userData.gapclosingFig = userData.fun_cost_gap{procID}('mainFig', handles.figure1, procID);
end
set(handles.figure1, 'UserData', userData);



function edit_maxgap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maxgap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maxgap as text
%        str2double(get(hObject,'String')) returns contents of edit_maxgap as a double


% --- Executes during object creation, after setting all properties.
function edit_maxgap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maxgap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_minlength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_minlength as text
%        str2double(get(hObject,'String')) returns contents of edit_minlength as a double


% --- Executes during object creation, after setting all properties.
function edit_minlength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_minlength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_merging.
function checkbox_merging_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_merging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_merging


% --- Executes on button press in checkbox_splitting.
function checkbox_splitting_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_splitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_splitting


% --- Executes on selection change in listbox_1.
function listbox_1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_1


% --- Executes during object creation, after setting all properties.
function listbox_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
% Hint: get(hObject,'Value') returns toggle state of checkbox_all
contents1 = get(handles.listbox_1, 'String');

chanIndex1 = get(handles.listbox_1, 'Userdata');
chanIndex2 = get(handles.listbox_2, 'Userdata');

% Return if listbox1 is empty
if isempty(contents1)
    return;
end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_2, 'String', contents1);
        chanIndex2 = chanIndex1;
    case 0
        set(handles.listbox_2, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
end
set(handles.listbox_2, 'UserData', chanIndex2);


% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents1 = get(handles.listbox_1, 'String');
contents2 = get(handles.listbox_2, 'String');
id = get(handles.listbox_1, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_1, 'Userdata');
chanIndex2 = get(handles.listbox_2, 'Userdata');

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_2, 'String', contents2, 'Userdata', chanIndex2);


% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Call back function of 'delete' button
contents = get(handles.listbox_2,'String');
id = get(handles.listbox_2,'Value');

% Return if list is empty
if isempty(contents) || isempty(id)
    return;
end

% Delete selected item
contents(id) = [ ];

% Delete userdata
chanIndex2 = get(handles.listbox_2, 'Userdata');
chanIndex2(id) = [ ];
set(handles.listbox_2, 'Userdata', chanIndex2);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_2,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_2,'String',contents);


% --- Executes on selection change in listbox_2.
function listbox_2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_2


% --- Executes during object creation, after setting all properties.
function listbox_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_verbose.
function checkbox_verbose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_verbose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_verbose



function edit_probDim_Callback(hObject, eventdata, handles)
% hObject    handle to edit_probDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_probDim as text
%        str2double(get(hObject,'String')) returns contents of edit_probDim as a double


% --- Executes during object creation, after setting all properties.
function edit_probDim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_probDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_histogram.
function checkbox_histogram_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_histogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_histogram


% --- Executes on button press in pushbutton_set_kalman.
function pushbutton_set_kalman_Callback(hObject, eventdata, handles)
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig - the handle of setting panel for kalman filter
%       initilization

userData = get(handles.figure1, 'UserData');
procID = get(handles.popupmenu_kalman_initialize, 'Value');
if procID > length(userData.fun_kalman_initialize)
    warndlg('Please select an option in the drop-down menu.','Error','modal')
    return
else
    userData.kalmanFig = userData.fun_kalman_initialize{procID}('mainFig', handles.figure1, procID);
end
set(handles.figure1, 'UserData', userData);


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

% Show the actual file name on GUI
i = userData.crtProc.channelIndex_(1);
file = ['Channel_' num2str(i) '_' file];
if ~userData.crtProc.overwrite_
    file = enumFileName(path, file);
end

str = [path file];
if length(str)>100
    str = ['...' str(end-100:end)];
end
set(handles.text_path, 'String', str)

set(handles.figure1, 'UserData', userData)


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


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig

userData = get(handles.figure1, 'UserData');

% Delete setting panel(single)
   
if isfield(userData, 'linkingFig') && ishandle(userData.linkingFig)
    
    delete(userData.linkingFig);
end

if isfield(userData, 'gapclosingFig') && ishandle(userData.gapclosingFig)
    
    delete(userData.gapclosingFig);
end

if isfield(userData, 'kalmanFig') && ishandle(userData.kalmanFig)
    
    delete(userData.kalmanFig);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
%       userData.linkingFig - the handle of setting panel for linking set-up
%       userData.gapclosingFig - the handle of setting panel for gap closing set-up
%       userData.kalmanFig

userData = get(handles.figure1, 'UserData');

fig(1) = isfield(userData, 'linkingFig') && ishandle(userData.linkingFig);
fig(2) = isfield(userData, 'gapclosingFig') && ishandle(userData.gapclosingFig);
fig(3) = isfield(userData, 'kalmanFig') && ishandle(userData.kalmanFig);

if any(fig)
    
    temp = find(fig);
    
    switch temp(1)
        
        case 1
            hPopupmenu = handles.popupmenu_linking;
        case 2
            hPopupmenu = handles.popupmenu_gapclosing;
        case 3
            hPopupmenu = handles.popupmenu_kalman_initialize;
    end
    str = get(hPopupmenu, 'String');
    user_response = questdlg(['The setting panel of method "',str{get(hPopupmenu, 'Value')},'" is still open. Do you want to close the setting panel of tracking?'], ...
        'Setting Open', 'No','Yes','Yes');
    
    if strcmpi(user_response, 'yes')

        delete(handles.figure1)
    end   
    
else

    delete(handles.figure1);
end


% --- Executes on button press in checkbox_overwrite.
function checkbox_overwrite_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_overwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_overwrite
