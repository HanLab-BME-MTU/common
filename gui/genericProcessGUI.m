function varargout = genericProcessGUI(varargin)
% genericProcessGUI M-file for genericProcessGUI.fig
%      genericProcessGUI, by itself, creates a new genericProcessGUI or raises the existing
%      singleton*.
%
%      H = genericProcessGUI returns the handle to a new genericProcessGUI or the handle to
%      the existing singleton*.
%
%      genericProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in genericProcessGUI.M with the given input arguments.
%
%      genericProcessGUI('Property','Value',...) creates a new genericProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before genericProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to genericProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help genericProcessGUI

% Last Modified by GUIDE v2.5 12-May-2011 17:12:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @genericProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @genericProcessGUI_OutputFcn, ...
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


% --- Executes just before genericProcessGUI is made visible.
function genericProcessGUI_OpeningFcn(hObject, ~, handles, varargin)
% Available tools 
% UserData data:
%       userData.mainFig - handle of main figure
%       userData.handles_main - 'handles' of main figure
%       userData.procID - The ID of process in the current package
%       userData.crtProc - handle of current process
%       userData.crtPackage - handles of current package
%       userData.procConstr - constructor of current process
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for genericProcessGUI
handles.output = hObject;

% Get main figure handle and process id
t = find(strcmp(varargin,'mainFig'));
userData.mainFig = varargin{t+1};
userData.procID = varargin{t+2};
userData.handles_main = guidata(userData.mainFig);

% Get current package and process
userData_main = get(userData.mainFig, 'UserData');
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get current process constructor
crtProcName = userData.crtPackage.processClassNames_{userData.procID};
userData.procConstr = str2func(crtProcName);
procString = [' Step ' num2str(userData.procID) ':' regexprep(crtProcName,'([A-Z])',' $1')];
set(handles.text_processName,'String',procString)
    
% If process does not exist, create a default one in user data.
if isempty(userData.crtProc)
    userData.crtProc = userData.procConstr(userData_main.MD(userData_main.id), ...
                                userData.crtPackage.outputDirectory_);
end

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% ---------------------- Channel Setup -------------------------

funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableChannels, 'String', {userData_main.MD(userData_main.id).channels_.channelPath_},...
        'Userdata', 1: length(userData_main.MD(userData_main.id).channels_));
    
% Set up selected input data channels and channel index
parentI = find( userData.crtPackage.depMatrix_(userData.procID,:) );

if isempty(parentI) || ~isempty( userData.crtPackage.processes_{userData.procID} )
    
    % If process has no dependency, or process already exists, display saved channels 
    set(handles.listbox_selectedChannels, 'String', ...
        {userData_main.MD(userData_main.id).channels_(funParams.ChannelIndex).channelPath_}, ...
        'Userdata',funParams.ChannelIndex);
    
elseif isempty( userData.crtPackage.processes_{userData.procID} )
    % If new process
        empty = false;
        for i = parentI
           if isempty(userData.crtPackage.processes_{i})
               empty = true;
               break;
           end
        end
            
        if ~empty

            % If all dependent processes exist
            channelIndex = userData.crtPackage.processes_{parentI(1)}.funParams_.ChannelIndex;
            for i = 2: length(parentI)
                channelIndex = intersect(channelIndex, ...
                    userData.crtPackage.processes_{parentI(i)}.funParams_.ChannelIndex);
            end  
            
            if ~isempty(channelIndex)
                set(handles.listbox_selectedChannels, 'String', ...
                    {userData_main.MD(userData_main.id).channels(channelIndex).channelPath_}, ...
                    'Userdata',channelIndex);    
            end
        end
end

% ---------------------- Parameter Setup -------------------------

% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
set(handles.figure1,'CurrentAxes',handles.axes_help);
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
function varargout = genericProcessGUI_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
% Delete figure
delete(handles.figure1);

% --- Executes on button press in checkbox_all.
function checkbox_all_Callback(hObject, ~, handles)

% Hint: get(hObject,'Value') returns toggle state of checkbox_all
availableChannels = get(handles.listbox_availableChannels, 'String');

availableChannelsIndx = get(handles.listbox_availableChannels, 'Userdata');
selectedChannelsIndx = get(handles.listbox_selectedChannels, 'Userdata');

% Return if listbox1 is empty
if isempty(availableChannels), return; end

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedChannels, 'String', availableChannels);
        selectedChannelsIndx = availableChannelsIndx;
    case 0
        set(handles.listbox_selectedChannels, 'String', {}, 'Value',1);
        selectedChannelsIndx = [ ];
end
set(handles.listbox_selectedChannels, 'UserData', selectedChannelsIndx);

% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button

availableChannels = get(handles.listbox_availableChannels, 'String');
selectedChannels = get(handles.listbox_selectedChannels, 'String');
id = get(handles.listbox_availableChannels, 'Value');

% If channel has already been added, return;
availableChannelsIndx = get(handles.listbox_availableChannels, 'Userdata');
selectedChannelsIndx = get(handles.listbox_selectedChannels, 'Userdata');

for i = id
    if any(strcmp(availableChannels{i}, selectedChannels) )
        continue;
    else
        selectedChannels{end+1} = availableChannels{i};
        selectedChannelsIndx = cat(2, selectedChannelsIndx, availableChannelsIndx(i));
    end
end

set(handles.listbox_selectedChannels, 'String', selectedChannels, 'Userdata', selectedChannelsIndx);

% --- Executes on button press in pushbutton_delete.
function pushbutton_delete_Callback(hObject, eventdata, handles)
% Call back function of 'delete' button
selectedChannels = get(handles.listbox_selectedChannels,'String');
id = get(handles.listbox_selectedChannels,'Value');

% Return if list is empty
if isempty(selectedChannels) || isempty(id),return; end

% Delete selected item
selectedChannels(id) = [ ];

% Delete userdata
selectedChannelsIndx = get(handles.listbox_selectedChannels, 'Userdata');
selectedChannelsIndx(id) = [ ];
set(handles.listbox_selectedChannels, 'Userdata', selectedChannelsIndx);

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(selectedChannels) && id>1)
    set(handles.listbox_selectedChannels,'Value',length(selectedChannels));
end
% Refresh listbox
set(handles.listbox_selectedChannels,'String',selectedChannels);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, ~, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

if isfield(userData, 'previewFig') && ishandle(userData.previewFig)
   delete(userData.previewFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(~, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
userData_main = get(userData.mainFig, 'UserData');

% -------- Check user input --------

if isempty(get(handles.listbox_selectedChannels, 'String'))
    errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal')
    return;
end

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
    threshold = get(handles.listbox_thresholdValues, 'String');
    if isempty(threshold)
       errordlg('Please provide at least one threshold value.','Setting Error','modal')
       return
    elseif length(threshold) ~= 1 && length(threshold) ~= length(get(handles.listbox_selectedChannels, 'String'))
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
   

% -------- Process Sanity check --------
% ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data.'],...
                'Setting Error','modal');
    return;
end

% -------- Set parameter --------

funParams = userData.crtProc.funParams_;
% Get parameter

% Set parameters
userData.crtProc.setPara(funParams);

% --------------------------------------------------

% If this is a brand new process, attach current process to MovieData and 
% package's process list 
if isempty( userData.crtPackage.processes_{userData.procID} )
    
    % Add new process to both process lists of MovieData and current package
    userData_main.MD(userData_main.id).addProcess( userData.crtProc );
    userData.crtPackage.setProcess(userData.procID, userData.crtProc);
    
    % Set font weight of process name bold
    set(userData.handles_main.(['checkbox_' num2str(userData.procID)]),...
            'FontWeight','bold');
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
        
        % Customize funParams to other movies
        % ChannelIndex - all channels
        % OutputDirectory - pacakge output directory
        
        l = length(userData_main.MD(x).channels_);
        temp = arrayfun(@(x)(x > l),channelIndex, 'UniformOutput', true );
        funParams.ChannelIndex = channelIndex(logical(~temp));
        
        if get(handles.checkbox_auto, 'value')
            
            funParams.ThresholdValue = [ ];
        else
            if length(threshold) == 1
                funParams.ThresholdValue = repmat(threshold, [1 length(funParams.ChannelIndex)]);
            else
                funParams.ThresholdValue = threshold(logical(~temp));
            end
        end
        
        funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'masks'];
        
        % if new process, create a new process with funParas and add to
        % MovieData and package's process list
        if isempty(userData_main.package(x).processes_{userData.procID})            
            process = userData.procConstr(userData_main.MD(x), userData_main.package(x).outputDirectory_, funParams);
            userData_main.MD(x).addProcess( process )
            userData_main.package(x).setProcess(userData.procID, process )            
        else
            % if process exist, replace the funParams with the new one
            userData_main.package(x).processes_{userData.procID}.setPara(funParams)
        end
        
        
        % Do sanity check - only check changed parameters
        procEx = userData_main.package(x).sanityCheck(false,'all');
        
        % Draw some bugs on the wall
        for i = 1: length(procEx)
            if ~isempty(procEx{i})
                % Record the icon and message to user data
                userData_main.statusM(x).IconType{i} = 'warn';
                userData_main.statusM(x).Msg{i} = procEx{i}(1).message;
            end
        end
    end
    
    % Save user data
    set(userData.mainFig, 'UserData', userData_main)

end
% -------------------------------------------------------------------------

% Save user data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);
