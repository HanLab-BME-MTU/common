function varargout = detectionProcessGUI(varargin)
% DETECTIONPROCESSGUI M-file for detectionProcessGUI.fig
%      DETECTIONPROCESSGUI, by itself, creates a new DETECTIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = DETECTIONPROCESSGUI returns the handle to a new DETECTIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      DETECTIONPROCESSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTIONPROCESSGUI.M with the given input arguments.
%
%      DETECTIONPROCESSGUI('Property','Value',...) creates a new DETECTIONPROCESSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detectionProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detectionProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detectionProcessGUI

% Last Modified by GUIDE v2.5 15-Dec-2010 15:12:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detectionProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @detectionProcessGUI_OutputFcn, ...
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


% --- Executes just before detectionProcessGUI is made visible.
function detectionProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.setFig(procID) = segmentProcessGUI('mainFig',handles.figure1, procID);
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
%
%       userData.segProc - cell array of segmentation processes, created
%                          after setting up segmentation processes
%
%
%       userData.procSetting - cell array of set-up GUIs of available
%                              processes
%       userData.procName - cell array of available segmentation processes
%       userData.procConstr - constructor of current process
%
%
%       userData.questIconData - help icon image information
%       userData.colormap - color map information
%
%       userData.set1Fig - the handle of setting panel for mask refinement
%                          process
%

[copyright openHelpFile] = userfcn_softwareConfig(handles);
set(handles.text_copyright, 'String', copyright)

userData = get(handles.figure1, 'UserData');
% Choose default command line output for segmentProcessGUI
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

% Get current process constructer, set-up GUIs and mask refinement process
% constructor
     
userData.procSetting = {@subResolutionDetectionGUI};
userData.procName = {'SubResolutionProcess'};                  
userData.procConstr = {@SubResolutionProcess};
popupMenuProcName = {'Sub-Resolution Object Detection',...
                     'Choose ...'};

% Initialize segProc and maskRefinProc in user data
userData.segProc = cell(1, length(userData.procName));

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% ---------------------- Channel and Parameter Setup  -------------------

set(handles.popupmenu_1, 'String', popupMenuProcName)

% Set up available input channels
set(handles.listbox_1, 'String', {userData.MD.channels_.channelPath_},...
        'Userdata', 1: length(userData.MD.channels_));
    
   
% Set up input channel list box
if isempty(userData.crtProc)
    
    set(handles.listbox_2, 'String', {userData.MD.channels_(1).channelPath_},...
        'Userdata', 1);
    
    % Set up pop-up menu
    set(handles.popupmenu_1, 'Value', length(get(handles.popupmenu_1, 'String')))
    
else
    
    channelIndex = userData.crtProc.channelIndex_;
    
    % If process has no dependency, or process already exists, display saved channels 
    set(handles.listbox_2, 'String', ...
        {userData.MD.channels_(channelIndex).channelPath_}, ...
        'Userdata',channelIndex);
    
    set(handles.popupmenu_1, 'Value', find(strcmp(userData.procName, class(userData.crtProc))) )
    set(handles.pushbutton_set_1, 'Enable', 'on')
%     set(handles.checkbox_overwrite, 'Value', userData.crtProc.overwrite_)
    
end

set(handles.checkbox_applytoall, 'Value', userData_main.applytoall(userData.procID));
 
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
    set(Img, 'UserData', struct('class', 'DetectionProcess'))
end



% ----------------------------------------------------------------

% Update user data and GUI data
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = detectionProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');
userData_main = get(userData.mainFig, 'UserData');

id = get(handles.popupmenu_1, 'Value');

if id == length(get(handles.popupmenu_1, 'String'))
    errordlg('Please select a method to segment your data.','Setting Error','modal')
    return 
end




sameprocess = false;

% Check if it needs to add a new segmentation process to movie data
if isempty(userData.crtProc) || ~isa(userData.crtProc, userData.procName{id})
    
    if isempty(userData.segProc{id})
        
        % Push user to set up the process
%         userData.crtProc = userData.procConstr{id}(userData.MD, userData.crtPackage.outputDirectory_);
%         userData.segProc{id} = userData.crtProc;
        
    else
        userData.crtProc = userData.segProc{id};
    end
else

    sameprocess = true;

end

% -------- Check user input --------

if isempty(get(handles.listbox_2, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
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

%---------Check if channel indexs are changed---------

channelIndex = get (handles.listbox_2, 'Userdata');
funParams = userData.crtProc.funParams_;

if sameprocess
    
    if ~isempty( setdiff(channelIndex, userData.crtProc.channelIndex_) ) ...
        || ~isempty( setdiff(userData.crtProc.channelIndex_, channelIndex) )

        % If channel indexs are changed, set procChanged to true
        userData.crtProc.setProcChanged(true);
    end
end
    
% Set channels
userData.crtProc.setChannelIndex(channelIndex)

% Remove background if settings are copied to other movies

if get(handles.checkbox_applytoall, 'Value')
    confirmApplytoAll = questdlg(['You are about to copy the current process settings to all movies.'...
        ' Previous settings as wellas  absolute background information will be lost. Do you want to continue?'],...
        'Apply settings to all movies','Yes','No','Yes');
    
    if strcmp(confirmApplytoAll,'Yes'),
        % Remove background information
        funParams = userData.crtProc.funParams_;
        funParams.detectionParam.background = [];
        userData.crtProc.setPara(funParams);
    else
        set(handles.checkbox_applytoall,'Value',0.0);            
        return
    end
end


% -------------------------- Assign process ----------------------------
% If this is not the original process, set the process to package

if ~sameprocess
    % Once changed method, delete original process in movie data's process
    % list
    if ~isempty( userData.crtPackage.processes_{userData.procID} )
        
        % Delete original process
        userData.MD.deleteProcess(userData.crtPackage.processes_{userData.procID})
        userData.crtPackage.setProcess(userData.procID, [ ])                                
    end
    
    % Add new process to both process lists of MovieData and current package
    userData.MD.addProcess( userData.crtProc );
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
        % 

        % Channel Index
        l = length(userData_main.MD(x).channels_);
        temp = arrayfun(@(x)(x > l),channelIndex, 'UniformOutput', true );
        channelIndex_ = channelIndex(logical(~temp));

        % Set output dir to default
        funParams.saveResults.dir = [userData_main.package(x).outputDirectory_  filesep userData.crtProc.name_];


        % if new process, create a new process with funParas and add to
        % MovieData and package's process list
        if isempty(userData_main.package(x).processes_{userData.procID})

            process = userData.procConstr{id}(userData_main.MD(x), userData_main.package(x).outputDirectory_, channelIndex_, funParams);
            process.setFileName(userData.crtProc.filename_)
            userData_main.MD(x).addProcess( process )
            userData_main.package(x).setProcess(userData.procID, process )

        % If process exists, same method, replace the funParams with the new one
        % If mask refinement exist, 
        elseif isa( userData_main.package(x).processes_{userData.procID}, userData.procName{id} )

            userData_main.package(x).processes_{userData.procID}.setPara(funParams)
            userData_main.package(x).processes_{userData.procID}.setChannelIndex(channelIndex_)
            userData_main.package(x).processes_{userData.procID}.setFileName(userData.crtProc.filename_)

        % if process exists, differenct method
        else

            % Delete segmentation process
            userData_main.MD(x).deleteProcess(userData_main.package(x).processes_{userData.procID})
            userData_main.package(x).setProcess(userData.procID, [ ])     

            % Add new segmentation process to package and movie data
            process = userData.procConstr{id}(userData_main.MD(x), userData_main.package(x).outputDirectory_, channelIndex_,funParams);
            process.setFileName(userData.crtProc.filename_)
            userData_main.MD(x).addProcess( process )
            userData_main.package(x).setProcess(userData.procID, process )       
       end


       % If current process is changed, then assume funParams are changed in
       % all movies
       if userData.crtProc.procChanged_ 

           userData_main.package(x).processes_{userData.procID}.setProcChanged(true);
       end

%            userData_main.package(x).processes_{userData.procID}.setOverwrite(get(handles.checkbox_overwrite, 'Value'))

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

end
% -------------------------------------------------------------------------

% Save user data
userData_main.applytoall(userData.procID)=get(handles.checkbox_applytoall,'Value');
set(userData.mainFig, 'UserData', userData_main)
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
delete(handles.figure1);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'set1Fig') && ishandle(userData.set1Fig)
    
    str = get(handles.popupmenu_1, 'String');
    user_response = questdlg(['The setting panel of method "',str{get(handles.popupmenu_1, 'Value')},'" is still open. Do you want to close the setting panel of detecion?'], ...
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


% --- Executes on selection change in popupmenu_1.
function popupmenu_1_Callback(hObject, eventdata, handles)

content = get(hObject, 'string');
if get(hObject, 'Value') == length(content)
    set(handles.pushbutton_set_1, 'Enable', 'off')
else
    set(handles.pushbutton_set_1, 'Enable', 'on')
end


% --- Executes during object creation, after setting all properties.
function popupmenu_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_set_1.
function pushbutton_set_1_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');
procID = get(handles.popupmenu_1, 'Value');
set1Fig = userData.procSetting{procID}('mainFig',handles.figure1,procID);
userData = get(handles.figure1, 'UserData');
userData.set1Fig = set1Fig;
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);



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
% call back function of 'select' button

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


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)


userData = get(handles.figure1, 'UserData');

% Delete setting panel(single)
   
if isfield(userData, 'set1Fig') && ishandle(userData.set1Fig)
    
    delete(userData.set1Fig);
    

    
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


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
userData = get(handles.figure1, 'UserData');

if isfield(userData, 'set1Fig') && ishandle(userData.set1Fig)
    
    str = get(handles.popupmenu_1, 'String');
    user_response = questdlg(['The setting panel of method "',str{get(handles.popupmenu_1, 'Value')},'" is still open. Do you want to close the setting panel of detecion?'], ...
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
