function varargout = thresholdProcessGUI(varargin)
% thresholdProcessGUI M-file for thresholdProcessGUI.fig
%      thresholdProcessGUI, by itself, creates a new thresholdProcessGUI or raises the existing
%      singleton*.
%
%      H = thresholdProcessGUI returns the handle to a new thresholdProcessGUI or the handle to
%      the existing singleton*.
%
%      thresholdProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in thresholdProcessGUI.M with the given input arguments.
%
%      thresholdProcessGUI('Property','Value',...) creates a new thresholdProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before thresholdProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to thresholdProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help thresholdProcessGUI

% Last Modified by GUIDE v2.5 01-Jun-2011 14:49:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @thresholdProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @thresholdProcessGUI_OutputFcn, ...
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


% --- Executes just before thresholdProcessGUI is made visible.
function thresholdProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:});


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
    
    if ~emptyParentProc
        parentChannelIndex = @(x) userData.crtPackage.processes_{x}.funParams_.ChannelIndex;
        channelIndex = 1:numel(userData.MD.channels_);
        for i = parentProc
            channelIndex = intersect(channelIndex,parentChannelIndex(i));
        end
        
        if ~isempty(channelIndex)
            channelString = userData.MD.getChannelPaths(channelIndex);
        else
            channelString = {};
        end
    end
end

set(handles.listbox_selectedChannels,'String',channelString,...
    'UserData',channelIndex);  

% ---------------------- Parameter Setup -------------------------

threshMethods = userData.crtProc.getMethods();
set(handles.popupmenu_thresholdingMethod,'String',{threshMethods(:).name},...
    'Value',funParams.MethodIndx);

if isempty(funParams.ThresholdValue)
   
    set(handles.checkbox_auto, 'Value', 1);
    set(get(handles.uipanel_automaticThresholding,'Children'),'Enable','on');
    set(get(handles.uipanel_fixedThreshold,'Children'),'Enable','off');
    if funParams.MaxJump
        set(handles.checkbox_max, 'Value', 1);
        set(handles.edit_jump, 'Enable','on','String',funParams.MaxJump);
    else
        set(handles.edit_jump, 'Enable', 'off');
    end
    set(handles.edit_GaussFilterSigma,'String',funParams.GaussFilterSigma);
    nSelectedChannels  = numel(get(handles.listbox_selectedChannels, 'String'));
    set(handles.listbox_thresholdValues, 'String',...
        num2cell(zeros(1,nSelectedChannels)));
    userData.thresholdValue=0;
else
    set(handles.checkbox_auto, 'Value', 0);
    set(get(handles.uipanel_fixedThreshold,'Children'),'Enable','on');
    set(get(handles.uipanel_automaticThresholding,'Children'),'Enable','off');
    set(handles.listbox_thresholdValues, 'String', num2cell(funParams.ThresholdValue));
    userData.thresholdValue=funParams.ThresholdValue(1);
end

% Save the image directories and names (for threshold preview)
userData.imageFileNames = userData.MD.getImageFileNames();
userData.imDirs  = userData.MD.getChannelPaths();

% Read the first image and update the sliders max value and steps
props=get(handles.listbox_selectedChannels,{'String','Value'});
userData.chanIndx = find(strcmp(props{1}{props{2}},userData.imDirs));
userData.imIndx=1;

% Initialize the frame number slider and eidt
nFrames=numel(userData.imageFileNames{userData.chanIndx});
set(handles.slider_frameNumber,'Value',userData.imIndx,'Min',1,...
    'Max',nFrames,'SliderStep',[1/double(nFrames)  10/double(nFrames)]);
set(handles.text_nFrames,'String',['/ ' num2str(nFrames)]);
set(handles.edit_frameNumber,'Value',userData.imIndx);

% Load the first image and update the threshold slide
userData.imData = imread([userData.imDirs{userData.chanIndx} filesep...
    userData.imageFileNames{userData.chanIndx}{userData.imIndx}]);
maxThresholdValue=max(max(userData.imData));
thresholdStep = 1/double(maxThresholdValue);

set(handles.edit_threshold,'String',userData.thresholdValue);
set(handles.slider_threshold,'Value',userData.thresholdValue,...
    'Max',maxThresholdValue,...
    'SliderStep',[thresholdStep  10*thresholdStep]);

% Update user data and GUI data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = thresholdProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% Delete figure
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');

% -------- Check user input --------

if isempty(get(handles.listbox_selectedChannels, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
end

if get(handles.checkbox_auto, 'value')
    if isnan(str2double(get(handles.edit_GaussFilterSigma, 'String'))) ...
            || str2double(get(handles.edit_GaussFilterSigma, 'String')) < 0
        errordlg(['Please provide a valid input for '''...
            get(handles.text_GaussFilterSigma,'String') '''.'],'Setting Error','modal');
        return;
    end
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
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

if get(handles.checkbox_auto, 'value')
    % if automatic thresholding
    funParams.ThresholdValue = [ ];
    funParams.MethodIndx=get(handles.popupmenu_thresholdingMethod,'Value');
    funParams.GaussFilterSigma = str2double(get(handles.edit_GaussFilterSigma,'String'));
    if get(handles.checkbox_max, 'value')
        funParams.MaxJump = str2double(get(handles.edit_jump,'String'));
    else
        funParams.MaxJump = 0;
    end
else
    funParams.ThresholdValue = threshold;
end

processGUI_ApplyFcn(hObject, eventdata, handles,funParams);

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

switch get(hObject,'Value')
    case 1
        set(handles.listbox_selectedChannels, 'String', contents1);
        chanIndex2 = chanIndex1;
        thresholdValues =zeros(1,numel(chanIndex1));
    case 0
        set(handles.listbox_selectedChannels, 'String', {}, 'Value',1);
        chanIndex2 = [ ];
        thresholdValues = [];
end
set(handles.listbox_selectedChannels, 'UserData', chanIndex2);
set(handles.listbox_thresholdValues,'String',num2cell(thresholdValues))
update_data(hObject,eventdata,handles);

% --- Executes on button press in pushbutton_select.
function pushbutton_select_Callback(hObject, eventdata, handles)
% call back function of 'select' button

contents1 = get(handles.listbox_availableChannels, 'String');
contents2 = get(handles.listbox_selectedChannels, 'String');
id = get(handles.listbox_availableChannels, 'Value');

% If channel has already been added, return;
chanIndex1 = get(handles.listbox_availableChannels, 'Userdata');
chanIndex2 = get(handles.listbox_selectedChannels, 'Userdata');
thresholdValues = cellfun(@str2num,get(handles.listbox_thresholdValues,'String'));

for i = id
    if any(strcmp(contents1{i}, contents2) )
        continue;
    else
        contents2{end+1} = contents1{i};
        thresholdValues(end+1) = 0;
        chanIndex2 = cat(2, chanIndex2, chanIndex1(i));

    end
end

set(handles.listbox_selectedChannels, 'String', contents2, 'Userdata', chanIndex2);
set(handles.listbox_thresholdValues,'String',num2cell(thresholdValues))
update_data(hObject,eventdata,handles);


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

% Update thresholdValues
thresholdValues = cellfun(@str2num,get(handles.listbox_thresholdValues,'String'));
thresholdValues(id) = [];
set(handles.listbox_thresholdValues,'String',num2cell(thresholdValues));

% Point 'Value' to the second last item in the list once the 
% last item has been deleted
if (id >length(contents) && id>1)
    set(handles.listbox_selectedChannels,'Value',length(contents));
    set(handles.listbox_thresholdValues,'Value',length(contents));
end
% Refresh listbox
set(handles.listbox_selectedChannels,'String',contents);
update_data(hObject,eventdata,handles);

% --- Executes on button press in checkbox_auto.
function checkbox_auto_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox_auto
switch get(hObject, 'Value')
    case 0
        set(get(handles.uipanel_automaticThresholding,'Children'),'Enable','off');
        set(get(handles.uipanel_fixedThreshold,'Children'),'Enable','on')
        set(handles.checkbox_max, 'Value', 0);
        set(handles.checkbox_applytoall, 'Value',0);
        update_data(hObject,eventdata,handles);
    case 1
        set(get(handles.uipanel_automaticThresholding,'Children'),'Enable','on');
        set(get(handles.uipanel_fixedThreshold,'Children'),'Enable','off'); 
        if ~get(handles.checkbox_max,'Value'), set(handles.edit_jump,'Enable','off'); end
        userData = get(handles.figure1, 'UserData');
        if ~isfield(userData, 'previewFig') || ishandle(userData.previewFig)
            delete(userData.previewFig);
        end
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
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

% --- Executes on button press in checkbox_max.
function checkbox_max_Callback(hObject, eventdata, handles)


switch get(hObject, 'value')
    case 0
        set(handles.edit_jump, 'Enable', 'off');
    case 1
        set(handles.edit_jump, 'Enable', 'on');
end


% --- Executes on key press with focus on pushbutton_done and none of its controls.
function pushbutton_done_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end

% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end


% --- Executes on button press in pushbutton_set_threshold.
function pushbutton_set_threshold_Callback(hObject, eventdata, handles)

newThresholdValue = get(handles.slider_threshold,'Value');
indx = get(handles.listbox_selectedChannels,'Value');
thresholdValues = cellfun(@str2num,get(handles.listbox_thresholdValues,'String'));
thresholdValues(indx) = newThresholdValue;
set(handles.listbox_thresholdValues,'String',num2cell(thresholdValues));


% --- Executes on button press in checkbox_preview.
function checkbox_preview_Callback(hObject, eventdata, handles)

if get(handles.checkbox_preview,'Value'), 
    update_data(hObject,eventdata,handles); 
else
    userData = get(handles.figure1, 'UserData');
    if isfield(userData, 'previewFig') && ishandle(userData.previewFig)
        delete(userData.previewFig);
    end
    % Save data and update graphics
    set(handles.figure1,'UserData',userData);
    guidata(hObject, handles);
end

function threshold_edition(hObject,eventdata, handles)

userData = get(handles.figure1, 'UserData');

% Retrieve the valuye of the selected threshold
if strcmp(get(hObject,'Tag'),'edit_threshold')
    thresholdValue = str2double(get(handles.edit_threshold, 'String'));
else
    thresholdValue = get(handles.slider_threshold, 'Value');
end
thresholdValue=round(thresholdValue);

% Check the validity of the supplied threshold
if isnan(thresholdValue) || thresholdValue < 0 || thresholdValue > get(handles.slider_threshold,'Max')
    warndlg('Please provide a valid coefficient.','Setting Error','modal');
    thresholdValue=userData.thresholdValue;
end

set(handles.slider_threshold,'Value',thresholdValue);
set(handles.edit_threshold,'String',thresholdValue);
    
% Save data and update graphics
set(handles.figure1,'UserData',userData);
guidata(hObject, handles);
update_data(hObject,eventdata,handles);


function imageNumber_edition(hObject,eventdata, handles)

% Retrieve the value of the selected image
if strcmp(get(hObject,'Tag'),'edit_imageNumber')
    imageNumber = str2double(get(handles.edit_frameNumber, 'String'));
else
    imageNumber = get(handles.slider_frameNumber, 'Value');
end
imageNumber=round(imageNumber);

% Check the validity of the supplied threshold
if isnan(imageNumber) || imageNumber < 0 || imageNumber > get(handles.slider_frameNumber,'Max')
    warndlg('Please provide a valid coefficient.','Setting Error','modal');
end

set(handles.slider_frameNumber,'Value',imageNumber);
set(handles.edit_frameNumber,'String',imageNumber);

% Save data and update graphics
guidata(hObject, handles);
update_data(hObject,eventdata,handles);


function update_data(hObject,eventdata, handles)

userData = get(handles.figure1, 'UserData');

if strcmp(get(get(hObject,'Parent'),'Tag'),'uipanel_2') || strcmp(get(hObject,'Tag'),'listbox_thresholdValues')
    % Check if changes have been at the list box level
    linkedListBoxes = {'listbox_2','listbox_thresholdValues'};
    checkLinkBox = strcmp(get(hObject,'Tag'),linkedListBoxes);
    if any(checkLinkBox)
        value = get(handles.(linkedListBoxes{checkLinkBox}),'Value');
        set(handles.(linkedListBoxes{~checkLinkBox}),'Value',value);
    else
        value = get(handles.listbox_selectedChannels,'Value');
    end
    thresholdString = get(handles.listbox_thresholdValues,'String');
    if ~isempty(thresholdString)
        thresholdValue = str2double(thresholdString{value});
    else
        thresholdValue=0;
    end
    set(handles.edit_threshold,'String',thresholdValue);
    set(handles.slider_threshold,'Value',thresholdValue);
    if isempty(thresholdString),
        if isfield(userData, 'previewFig'), delete(userData.previewFig); end
        set(handles.figure1, 'UserData', userData);
        guidata(hObject,handles);
        return
    end
end


% Retrieve the channex index
props=get(handles.listbox_selectedChannels,{'String','Value'});
chanIndx = find(strcmp(props{1}{props{2}},userData.imDirs));
imIndx = get(handles.slider_frameNumber,'Value');
thresholdValue = get(handles.slider_threshold, 'Value');

% Load a new image in case the image number or channel has been changed
if (chanIndx~=userData.chanIndx) ||  (imIndx~=userData.imIndx)
    userData.imData=imread([userData.imDirs{chanIndx} filesep...
        userData.imageFileNames{chanIndx}{imIndx}]);
    
    % Get the value of the new maximum threshold
    maxThresholdValue=max(max(userData.imData));
    % Update the threshold Value if above th new maximum
    thresholdValue=min(thresholdValue,maxThresholdValue);
    
    set(handles.slider_threshold,'Value',thresholdValue,'Max',maxThresholdValue,...
        'SliderStep',[1/double(maxThresholdValue)  10/double(maxThresholdValue)]);
    set(handles.edit_threshold,'String',thresholdValue);
    
    userData.updateImage=1;
    userData.chanIndx=chanIndx;
    userData.imIndx=imIndx;
else
    userData.updateImage=0;
end

% Check if the threshold map should be refreshed
if  (thresholdValue~=userData.thresholdValue);
    userData.thresholdValue = thresholdValue;
    userData.updateThreshold=1;
else
    userData.updateThreshold=0;
end

% Save the data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);

% Update graphics if applicable
if get(handles.checkbox_preview,'Value') && ~get(handles.checkbox_auto,'Value'),

    % Create figure if non-existing or closed
    if ~isfield(userData, 'previewFig') || ~ishandle(userData.previewFig)
        userData.previewFig = figure;
        userData.newFigure = 1;
    else
        figure(userData.previewFig);
        userData.newFigure = 1;
    end
    
    % Retrieve the handle of the figures image
    imHandle =get(get(userData.previewFig,'Children'),'Children');
    if userData.newFigure || userData.updateImage || isempty(imHandle)
        imagesc(userData.imData);
        axis off;
        imHandle =get(get(userData.previewFig,'Children'),'Children');
    end
    
    % Preview the tresholding output using the alphaData mapping
    if userData.newFigure || userData.updateThreshold
        alphamask=ones(size(userData.imData));
        alphamask(userData.imData<=userData.thresholdValue)=.4;
        set(imHandle,'AlphaData',alphamask,'AlphaDataMapping','none');
    end
    
    set(handles.figure1, 'UserData', userData);
    guidata(hObject,handles);
end
