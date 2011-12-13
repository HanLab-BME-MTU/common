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

% Last Modified by GUIDE v2.5 13-Dec-2011 17:58:33

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

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',1);

% Parameter Setup
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

set(handles.popupmenu_probDim,'String',{'2','3'},'UserData',[2 3],...
    'Value',find(funParams.probDim==[2 3]));
set(handles.checkbox_verbose, 'Value', funParams.verbose)

% gapCloseParam
set(handles.edit_maxgap, 'String', num2str(funParams.gapCloseParam.timeWindow - 1))
set(handles.edit_minlength, 'String', num2str(funParams.gapCloseParam.minTrackLen))
set(handles.checkbox_histogram, 'Value', funParams.gapCloseParam.diagnostics)

set(handles.checkbox_merging, 'Value',ismember(funParams.gapCloseParam.mergeSplit,[1 2]));
set(handles.checkbox_splitting, 'Value',ismember(funParams.gapCloseParam.mergeSplit,[1 3]));
    
% Set cost matrics
defaultLinkingCostMat = TrackingProcess.getDefaultLinkingCostMatrices(userData.MD,5);
defaultGapClosingCostMat = TrackingProcess.getDefaultGapClosingCostMatrices(userData.MD,5);
userData.cost_linking = {defaultLinkingCostMat.funcName};
userData.cost_gapclosing = {defaultGapClosingCostMat.funcName};
userData.fun_cost_linking = {defaultLinkingCostMat.GUI};
userData.fun_cost_gap = {defaultGapClosingCostMat.GUI};

% Retrieve index of default cost matrices
i1 = find(strcmp(funParams.costMatrices(1).funcName, userData.cost_linking));
i2 = find(strcmp(funParams.costMatrices(2).funcName, userData.cost_gapclosing));
assert(isscalar(i1) && isscalar(i2),'User-defined: the length of matching methods must be 1.')
u1 = cell(1, numel(defaultLinkingCostMat));
u2 = cell(1,numel(defaultGapClosingCostMat));
u1{i1} = funParams.costMatrices(1).parameters;
u2{i2} = funParams.costMatrices(2).parameters;

set(handles.popupmenu_linking, 'Value', i1, 'UserData', u1,...
    'String',{defaultLinkingCostMat.name})
set(handles.popupmenu_gapclosing, 'Value', i2, 'UserData', u2,...
    'String',{defaultGapClosingCostMat.name})


% Kalman functions
userData.reserveMemFunctions = TrackingProcess.getKalmanReserveMemFunctions;
userData.initializeFunctions = TrackingProcess.getKalmanInitializeFunctions;
userData.calcGainFunctions = TrackingProcess.getKalmanCalcGainFunctions;
userData.timeReverseFunctions = TrackingProcess.getKalmanTimeReverseFunctions;

i1 = find(strcmp(funParams.kalmanFunctions.reserveMem, {userData.reserveMemFunctions.funcName}));
i2 = find(strcmp(funParams.kalmanFunctions.initialize, {userData.initializeFunctions.funcName}));
i3 = find(strcmp(funParams.kalmanFunctions.calcGain, {userData.calcGainFunctions.funcName}));
i4 = find(strcmp(funParams.kalmanFunctions.timeReverse, {userData.timeReverseFunctions.funcName}));

assert(isscalar(i1) && isscalar(i2) && isscalar(i3) && isscalar(i4),...
    'User-defined: the length of matching methods must be 1.');

u2 = cell(1, numel(userData.initializeFunctions));
u2{i2} = funParams.costMatrices(1).parameters.kalmanInitParam;

set(handles.popupmenu_kalman_reserve, 'String', {userData.reserveMemFunctions.name}, 'Value', i1)
set(handles.popupmenu_kalman_initialize,'String', {userData.initializeFunctions.name}, 'Value', i2, 'UserData', u2)
set(handles.popupmenu_kalman_gain, 'String', {userData.calcGainFunctions.name}, 'Value', i3)
set(handles.popupmenu_kalman_reverse,'String', {userData.timeReverseFunctions.name}, 'Value', i4)

set(handles.checkbox_export, 'Value', funParams.saveResults.export)

% funParams.saveResults
props=get(handles.listbox_selectedChannels,{'UserData','Value'});
chanIndx=props{1}(props{2});
fullPath= [funParams.OutputDirectory filesep...
    'Channel_' num2str(chanIndx) '_' funParams.OutputFilename];
set(handles.text_path, 'String', fullPath)

% Initialize children figure handles
userData.linkingFig=-1;
userData.gapclosingFig=-1;
userData.kalmanFig=-1;

% Choose default command line output for trackingProcessGUI
handles.output = hObject;

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

% Check User Input
if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

props = get(handles.popupmenu_probDim, {'UserData','Value'});
probDim=props{1}(props{2});

timeWindow = str2double(get(handles.edit_maxgap, 'String'))+1;
if isnan(timeWindow) || timeWindow < 0 || floor(timeWindow) ~= ceil(timeWindow)
    errordlg('Please provide a valid value to parameter "Maximum Gap to Close".','Error','modal')
    return
end

minTrackLen = str2double(get(handles.edit_minlength, 'String'));
if isnan(minTrackLen) || minTrackLen < 0 || floor(minTrackLen) ~= ceil(minTrackLen)
    errordlg('Please provide a valid value to parameter "Minimum Length of Track Segment from First Step to use in Second Step".','Error','modal')
    return
end

% -------- Set parameter --------
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

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

funParams.saveResults.export = get(handles.checkbox_export, 'Value');

% Cost matrices
i_linking = get(handles.popupmenu_linking, 'Value');
i_gapclosing = get(handles.popupmenu_gapclosing, 'Value');

u_linking = get(handles.popupmenu_linking, 'UserData');
u_gapclosing = get(handles.popupmenu_gapclosing, 'UserData');

if isempty( u_linking{i_linking} )
    errordlg('Plese set up the selected cost function for "Step 1: frame-to-frame linking".','Error','modal')
end

if isempty( u_gapclosing{i_gapclosing} )
    errordlg('Plese set up the selected cost function for "Step 2: gap closing, mergin and splitting".','Error','modal')
end

funParams.costMatrices(1).funcName = userData.cost_linking{i_linking};
funParams.costMatrices(1).parameters = u_linking{i_linking};
funParams.costMatrices(2).funcName = userData.cost_gapclosing{i_gapclosing};
funParams.costMatrices(2).parameters = u_gapclosing{i_gapclosing};

% Get Kalman values
props = get(handles.popupmenu_kalman_initialize, {'Value','UserData'});
funParams.kalmanFunctions.initialize = userData.initializeFunctions(props{1}).funcName;
funParams.costMatrices(1).parameters.kalmanInitParam = props{2}{props{1}};
i = get(handles.popupmenu_kalman_reserve, 'Value');
funParams.kalmanFunctions.reserveMem  = userData.reserveMemFunctions(i).funcName;
i = get(handles.popupmenu_kalman_gain, 'Value');
funParams.kalmanFunctions.calcGain    = userData.calcGainFunctions(i).funcName;
i = get(handles.popupmenu_kalman_reverse, 'Value');
funParams.kalmanFunctions.timeReverse = userData.timeReverseFunctions(i).funcName;

% Set up parameters effected by funParams.gapCloseParam.timeWindow
funParams.costMatrices(2).parameters.brownStdMult = funParams.costMatrices(2).parameters.brownStdMult(1) * ones(funParams.gapCloseParam.timeWindow,1);
funParams.costMatrices(2).parameters.linStdMult = funParams.costMatrices(2).parameters.linStdMult(1) * ones(funParams.gapCloseParam.timeWindow,1);

processGUI_ApplyFcn(hObject,eventdata,handles,funParams)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure1);

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

maxgap = str2double(get(handles.edit_maxgap, 'String'));
if isnan(maxgap) || maxgap < 0 || floor(maxgap) ~= ceil(maxgap)
    errordlg('Please provide a valid value to parameter "Maximum Gap to Close".','Warning','modal')

else
    timeWindow = maxgap + 1; % Retrieve the new value for the time window

    % Retrieve the parameters of the linking and gap closing matrices
    u_linking = get(handles.popupmenu_linking, 'UserData');
    linkingID = get(handles.popupmenu_linking, 'Value');
    linkingParameters = u_linking{linkingID};
    u_gapclosing = get(handles.popupmenu_gapclosing, 'UserData');
    gapclosingID = get(handles.popupmenu_gapclosing, 'Value');
    gapclosingParameters = u_gapclosing{gapclosingID};

    % Check for changes
    linkingnnWindowChange=(linkingParameters.nnWindow~=timeWindow);
    gapclosingnnWindowChange=(gapclosingParameters.nnWindow~=timeWindow);
    gapclosingtimeReachConfBChange=(gapclosingParameters.timeReachConfB~=timeWindow);
    gapclosingtimeReachConfLChange=(gapclosingParameters.timeReachConfL~=timeWindow);

    if linkingnnWindowChange || gapclosingnnWindowChange ||...
            gapclosingtimeReachConfBChange || gapclosingtimeReachConfLChange
        % Optional: asks the user if the time window value should be propagated
        % to the linking and gap closing matrics
        modifyParameters=questdlg('Do you want to propagate the changes in the maximum number of gaps to close?',...
           'Parameters update','Yes','No','Yes');
        if strcmp(modifyParameters,'Yes')
            % Save changes
            linkingParameters.nnWindow=timeWindow;
            gapclosingParameters.nnWindow=timeWindow;
            gapclosingParameters.timeReachConfB=timeWindow;
            gapclosingParameters.timeReachConfL=timeWindow;
            
            u_linking{linkingID} = linkingParameters;
            u_gapclosing{gapclosingID} = gapclosingParameters;
            
            set(handles.popupmenu_linking, 'UserData', u_linking)
            set(handles.popupmenu_gapclosing, 'UserData', u_gapclosing)
            guidata(hObject,handles);
        end
    end
end

% --- Executes on button press in pushbutton_set_kalman.
function pushbutton_set_kalman_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
funcId = get(handles.popupmenu_kalman_initialize, 'Value');

if funcId > numel(userData.initializeFunctions)
    warndlg('Please select an option in the drop-down menu.','Error','modal')
    return
else
    userData.kalmanFig = userData.initializeFunctions(funcId).GUI('mainFig', handles.figure1, funcId);
end
set(handles.figure1, 'UserData', userData);


% --- Executes on button press in pushbutton_saveas.
function pushbutton_saveas_Callback(hObject, eventdata, handles)


% Ask user where to save the movie data file
fullPath = get(handles.text_path, 'String');
[file path] = uiputfile('*.mat','Save Results As',fullPath);
        
if ~any([file path]), return; end

tokens = regexp(file,'^Channel_\d_(.+)$','tokens');
if isempty(tokens)
    props=get(handles.listbox_selectedChannels,{'UserData','Value'});
    chanIndx=props{1}(props{2});
    fullPath= [path  'Channel_' num2str(chanIndx) '_' file];
else
    fullPath=[path file];
end

set(handles.text_path, 'String', fullPath)

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

% Delete setting panels
if ishandle(userData.linkingFig), delete(userData.linkingFig);end
if ishandle(userData.gapclosingFig), delete(userData.gapclosingFig); end
if ishandle(userData.kalmanFig), delete(userData.kalmanFig); end


% --- Executes on button press in checkbox_export.
function checkbox_export_Callback(hObject, eventdata, handles)

if get(hObject,'Value')
    exportMsg=sprintf('The output matrices resulting from this process might be very large. Be cautious if you have large movies');
    if any([get(handles.checkbox_merging, 'Value') get(handles.checkbox_splitting, 'Value')])
        exportMsg =[exportMsg sprintf('\n \nAny merging and splitting information will be lost in the exported format.')];
    end
    warndlg(exportMsg,'Warning','modal')
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)

delete(hObject);
