function varargout = maskRefinementProcessGUI(varargin)
% MASKREFINEMENTPROCESSGUI M-file for maskRefinementProcessGUI.fig
%      MASKREFINEMENTPROCESSGUI, by itself, creates a new MASKREFINEMENTPROCESSGUI or raises the existing
%      singleton*.
%
%      H = MASKREFINEMENTPROCESSGUI returns the handle to a new MASKREFINEMENTPROCESSGUI or the handle to
%      the existing singleton*.
%
%      MASKREFINEMENTPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASKREFINEMENTPROCESSGUI.M with the given input arguments.
%
%      MASKREFINEMENTPROCESSGUI('Property','Value',...) creates a new MASKREFINEMENTPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before maskRefinementProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to maskRefinementProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help maskRefinementProcessGUI

% Last Modified by GUIDE v2.5 01-Jun-2011 16:59:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @maskRefinementProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @maskRefinementProcessGUI_OutputFcn, ...
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


% --- Executes just before maskRefinementProcessGUI is made visible.
function maskRefinementProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)


processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},...
    @MaskRefinementProcess);

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;


% Set up available input channels
set(handles.listbox_availableChannels,...
    'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));  
    
% Set up selected input data channels and channel index
parentProc = find(userData.crtPackage.depMatrix_(userData.procID,:));

if isempty(parentProc) || ~isempty(userData.crtPackage.processes_{userData.procID})
    % If process has no dependency, or process already exists
    channelString =userData.MD.getChannelPaths(funParams.ChannelIndex);
    channelIndex = funParams.ChannelIndex;
        
elseif isempty(userData.crtPackage.processes_{userData.procID})
    % Check existence of all parent processes
    emptyParentProc = any(cellfun(@isempty,...
        userData.crtPackage.processes_(parentProc)));
    channelIndex = 1:numel(userData.MD.channels_);
    if ~emptyParentProc
        parentChannelIndex = @(x) userData.crtPackage.processes_{x}.funParams_.ChannelIndex;
        for i = parentProc
            channelIndex = intersect(channelIndex,parentChannelIndex(i));
        end
    end
    if ~isempty(channelIndex)
        channelString = userData.MD.getChannelPaths(channelIndex);
    else
        channelString = {};
    end
end

set(handles.listbox_selectedChannels,'String',channelString,...
    'UserData',channelIndex);  

    
% ---------------------- Parameter Setup -----------------------

if funParams.MaskCleanUp
    if ~funParams.FillHoles
        set(handles.checkbox_fillholes, 'Value', 0)
    end
    set(handles.edit_1, 'String',num2str(funParams.MinimumSize))
    set(handles.edit_2, 'String',num2str(funParams.ClosureRadius))
    set(handles.edit_3, 'String',num2str(funParams.ObjectNumber))
else
    set([handles.checkbox_cleanup handles.checkbox_fillholes], 'Value', 0);
    set(get(handles.uipanel_cleanup,'Children'),'Enable','off');
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
    

% Update user data and GUI data
handles.output = hObject;
set(handles.figure1, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = maskRefinementProcessGUI_OutputFcn(hObject, eventdata, handles) 
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
delete(handles.figure1);


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

% -------- Process Sanity check --------
% ( only check underlying data )

try
    userData.crtProc.sanityCheck;
catch ME

    errordlg([ME.message 'Please double check your data'],...
                'Setting Error','modal');
    return;
end

% -------- Set parameter --------

channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams = userData.crtProc.funParams_;
funParams.ChannelIndex = channelIndex;

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


% --------------------------------------------------

% If this is a brand new process, attach current process to MovieData and 
% package's process list 
if isempty( userData.crtPackage.processes_{userData.procID} )
    
    % Add new process to both process lists of MovieData and current package
    userData_main.MD(userData_main.id).addProcess( userData.crtProc );
    userData.crtPackage.setProcess(userData.procID, userData.crtProc);
    
    % Customized to Biosensors Package (hard coded)
    userData.crtPackage.setDepMatrix(7, userData.procID, 1);
    userData.crtPackage.setDepMatrix(9, userData.procID, 1);        
    
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
       userfcn_drawIcon(userData.handles_main,'warn',i,procEx{i}(1).message, true); % user data is retrieved, updated and submitted
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
   
   funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep 'refined_masks'];
   
   % if new process, create a new process with funParas and add to
   % MovieData and package's process list
   if isempty(userData_main.package(x).processes_{userData.procID})
       
       process = userData.procConstr(userData_main.MD(x), userData_main.package(x).outputDirectory_, funParams);
       userData_main.MD(x).addProcess( process )
       userData_main.package(x).setProcess(userData.procID, process )
       
   % if process exist, replace the funParams with the new one
   else
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

switch get(handles.checkbox_all,'Value')
    case 1
        set(handles.listbox_selectedChannels, 'String', contents1);
        chanIndex2 = chanIndex1;
    case 0
        set(handles.listbox_selectedChannels, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
end
set(handles.listbox_selectedChannels, 'UserData', chanIndex2);


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


% --- Executes on button press in checkbox_cleanup.
function checkbox_cleanup_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_auto

if get(hObject, 'Value')
    set(get(handles.uipanel_cleanup,'Children'),'Enable','on');
else
    set(get(handles.uipanel_cleanup,'Children'),'Enable','off');
end


% --- Executes on button press in checkbox_edge.
function checkbox_edge_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_auto

if get(hObject, 'Value')
    set(get(handles.uipanel_edge,'Children'),'Enable','on');
else
    set(get(handles.uipanel_edge,'Children'),'Enable','off');
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'helpFig') && ishandle(userData.helpFig)
   delete(userData.helpFig) 
end

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


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


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
