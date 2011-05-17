function varargout = segmentationProcessGUI(varargin)
%SEGMENTATIONPROCESSGUI M-file for segmentationProcessGUI.fig
%      SEGMENTATIONPROCESSGUI, by itself, creates a new SEGMENTATIONPROCESSGUI or raises the existing
%      singleton*.
%
%      H = SEGMENTATIONPROCESSGUI returns the handle to a new SEGMENTATIONPROCESSGUI or the handle to
%      the existing singleton*.
%
%      SEGMENTATIONPROCESSGUI('Property','Value',...) creates a new SEGMENTATIONPROCESSGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to segmentationProcessGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SEGMENTATIONPROCESSGUI('CALLBACK') and SEGMENTATIONPROCESSGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SEGMENTATIONPROCESSGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segmentationProcessGUI

% Last Modified by GUIDE v2.5 19-Apr-2011 16:54:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segmentationProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @segmentationProcessGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before segmentationProcessGUI is made visible.
function segmentationProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% userData.setFig(procID) = segmentationProcessGUI('mainFig',handles.figure1, procID);
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
%       userData.set1Fig - the handle of setting panel for segmentation
%                          process
%       userData.set2Fig - the handle of setting panel for mask refinement
%                          process
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
userData.MD = userData_main.MD(userData_main.id);  % Get the current Movie Data
userData.crtPackage = userData_main.crtPackage;
userData.crtProc = userData.crtPackage.processes_{userData.procID};

% Get current process constructor
crtProcName = userData.crtPackage.processClassNames_{userData.procID};
userData.procConstr = str2func(crtProcName);
procString = [' Step ' num2str(userData.procID) ':' regexprep(crtProcName,'([A-Z])',' $1')];
set(handles.text_processName,'String',procString);
figString = [' Setting - ' regexprep(crtProcName,'([A-Z])',' $1')];
set(handles.figure1,'Name',figString);

% Get current process constructer, set-up GUIs and mask refinement process
% constructor
     

userData.procSetting = {@thresholdProcessGUI};
userData.procName = {'ThresholdProcess'};                  
userData.procConstr = {@ThresholdProcess};
popupMenuProcName = {'Thresholding Segmentation',...
                     'Choose ...'};

% Initialize segProc in user data
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
    
    funParams = userData.crtProc.funParams_;
    
    % If process has no dependency, or process already exists, display saved channels 
    set(handles.listbox_2, 'String', ...
        {userData.MD.channels_(funParams.ChannelIndex).channelPath_}, ...
        'Userdata',funParams.ChannelIndex);
    

    set(handles.popupmenu_1, 'Value', find(strcmp(userData.procName, class(userData.crtProc))) )

    
    set(handles.pushbutton_set_1, 'Enable', 'on')
    
    % Set up post-processing
    if isempty(userData.crtProc.maskRefineProcess_)
        
       set(handles.checkbox_post, 'Value', 0)
       set(handles.pushbutton_set_2, 'Enable', 'off')
    end
    
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
    set(Img, 'UserData', struct('class',class(userData.crtPackage)))
end



% ----------------------------------------------------------------

% Update user data and GUI data
set(hObject, 'UserData', userData);

uicontrol(handles.pushbutton_done);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = segmentationProcessGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_preview.
function pushbutton_preview_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
        userData.crtProc = userData.procConstr{id}(userData.MD, userData.crtPackage.outputDirectory_);
        userData.segProc{id} = userData.crtProc;
        
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

channelIndex = get (handles.listbox_2, 'Userdata');

if ~isempty(userData.crtProc.funParams_.ThresholdValue) && length( userData.crtProc.funParams_.ChannelIndex ) ~= length(channelIndex)
    errordlg('The number of threshold is inconsistent with the number of channels.','Setting Error','modal') 
   return 
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
funParams = userData.crtProc.funParams_;
threshold = funParams.ThresholdValue;
   
% Get parameter
funParams.ChannelIndex = channelIndex;

% Set parameters
userData.crtProc.setPara(funParams);

% --------------------------------------------------


% If this is not the original process, set the process to package

if ~sameprocess
    % Once changed method, delete original process in movie data's process
    % list
    if ~isempty( userData.crtPackage.processes_{userData.procID} )
        
        % Delete mask refinement process attached to the segmentation
        % process, if there is any
        if ~isempty(userData.crtPackage.processes_{userData.procID}.maskRefineProcess_)
            
            userData.MD.deleteProcess(userData.crtPackage.processes_{userData.procID}.maskRefineProcess_)
            userData.crtPackage.processes_{userData.procID}.clearMaskRefineProcess 
        end
        
        % Delete segmentation process
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

% -------------- Attach Mask Refinement to Process -----------------------

% If do post processing
if get(handles.checkbox_post, 'Value')
    
    % If no maskrefinement process attached
    if isempty(userData.crtProc.maskRefineProcess_)
        
        if isempty(userData.maskRefineProc)
            % Create default mask refinement process and attach to
            % segmentation process
            userData.crtProc.setMaskRefineProcess( ...
                        userData.maskRefineConstr(userData.MD, userData.crtPackage.outputDirectory_)  )
            
        else
            % Link mask refinement process object in user data to
            % segmentation process
            
            userData.crtProc.setMaskRefineProcess( userData.maskRefineProc )
            
        end
        
        % Add mask refinement process to process lists of MovieData 
        userData.MD.addProcess( userData.crtProc.maskRefineProcess_ )
        
    end
    
    funParams_MR = userData.crtProc.maskRefineProcess_.funParams_;
    
    if ~isempty(userData.maskRefineProc)
       funParams_MR = userData.maskRefineProc.funParams_;
    end

    % Set Channel Index for mask refinement process
    
    funParams_MR.ChannelIndex = channelIndex;
    
    userData.crtProc.maskRefineProcess_.setPara(funParams_MR);
    
    
    
% If do not do post processing    
elseif ~isempty(userData.crtProc.maskRefineProcess_)
    
    % Delete mask refinement process attached
    userData.MD.deleteProcess(userData.crtProc.maskRefineProcess_)
    userData.crtProc.clearMaskRefineProcess
        
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
   % OutputDirectory - pacakge output directory
       
   l = length(userData_main.MD(x).channels_);
   temp = arrayfun(@(x)(x > l),channelIndex, 'UniformOutput', true );
   funParams.ChannelIndex = channelIndex(logical(~temp));
   
   if ~isempty(funParams.ThresholdValue)
       
        if length(threshold) == 1
            funParams.ThresholdValue = repmat(threshold, [1 length(funParams.ChannelIndex)]);
        else
            funParams.ThresholdValue = threshold(logical(~temp));
        end
   end
   
   funParams.OutputDirectory  = [userData_main.package(x).outputDirectory_  filesep userData.crtProc.name_];
   
   % Set mask refinement's parameter
   % ChannelIndex - all channels
   % OutputDirectory - pacakge output directory   
   if get(handles.checkbox_post, 'Value')
       
       funParams_MR.ChannelIndex = channelIndex(logical(~temp));
       funParams_MR.OutputDirectory = [userData_main.package(x).outputDirectory_  filesep 'refined_masks'];
   end
   
   % if new process, create a new process with funParas and add to
   % MovieData and package's process list
   if isempty(userData_main.package(x).processes_{userData.procID})
       
       process = userData.procConstr{id}(userData_main.MD(x), userData_main.package(x).outputDirectory_, funParams);
       userData_main.MD(x).addProcess( process )
       userData_main.package(x).setProcess(userData.procID, process )
       
   % If process exists, same method, replace the funParams with the new one
   % If mask refinement exist, 
   elseif isa( userData_main.package(x).processes_{userData.procID}, userData.procName{id} )
       userData_main.package(x).processes_{userData.procID}.setPara(funParams)
       
   % if process exists, differenct method
   else
       
        % Delete mask refinement process attached to the segmentation
        % process, if there is any
        if ~isempty(userData_main.package(x).processes{userData.procID}.maskRefineProcess_)
            
            userData_main.MD(x).deleteProcess(userData_main.package(x).processes_{userData.procID}.maskRefineProcess_)
            userData_main.package(x).processes_{userData.procID}.clearMaskRefineProcess            
        end       
       
        % Delete segmentation process
        userData_main.MD(x).deleteProcess(userData_main.package(x).processes_{userData.procID})
        userData_main.package(x).setProcess(userData.procID, [ ])     
       
        % Add new segmentation process to package and movie data
        process = userData.procConstr{id}(userData_main.MD(x), userData_main.package(x).outputDirectory_, funParams);
        userData_main.MD(x).addProcess( process )
        userData_main.package(x).setProcess(userData.procID, process )       
   end
   
   
   % Attach mask refinement process to segmentation process
   if get(handles.checkbox_post, 'Value')
       
       if isempty(userData_main.package(x).processes_{userData.procID}.maskRefineProcess_)
           % If no mask refinement process attached, create new mask refinement process and add to segmentation process
           % and movie data process list

           userData_main.package(x).processes_{userData.procID}.setMaskRefineProcess( ...
               userData.maskRefineConstr(userData_main.MD(x), userData_main.package(x).outputDirectory_, funParams_MR) )
           userData_main.MD(x).addProcess(userData_main.package(x).processes_{userData.procID}.maskRefineProcess_)

       else
           userData_main.package(x).processes_{userData.procID}.maskRefineProcess_.setPara(funParams_MR)
       end       
       
   % If no post processing selected, check if segmentation process has mask
   % refinement process attached to it
   
   elseif ~isempty(userData_main.package(x).processes_{userData.procID}.maskRefineProcess_)
       
   % Delete and clear attached mask refinement process    
           userData_main.MD(x).deleteProcess(userData_main.package(x).processes_{userData.procID}.maskRefineProcess_)
           userData_main.package(x).processes_{userData.procID}.clearMaskRefineProcess  
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


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
delete(handles.figure1);


% --- Executes on selection change in popupmenu_1.
function popupmenu_1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_1
content = get(hObject, 'string');
if get(hObject, 'Value') == length(content)
    set(handles.pushbutton_set_1, 'Enable', 'off')
else
    set(handles.pushbutton_set_1, 'Enable', 'on')
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

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');

% Delete setting panel(single)
if isfield(userData, 'set1Fig') && ishandle(userData.set1Fig)
   delete(userData.set1Fig) 
end

% Delete post-processing panel(single)
if isfield(userData, 'set2Fig') && ishandle(userData.set2Fig)
   delete(userData.set2Fig) 
end
