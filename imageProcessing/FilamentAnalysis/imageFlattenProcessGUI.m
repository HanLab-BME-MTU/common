function varargout = imageFlattenProcessGUI(varargin)
% imageFlattenProcessGUI M-file for imageFlattenProcessGUI.fig
%      imageFlattenProcessGUI, by itself, creates a new imageFlattenProcessGUI or raises the existing
%      singleton*.
%
%      H = imageFlattenProcessGUI returns the handle to a new imageFlattenProcessGUI or the handle to
%      the existing singleton*.
%
%      imageFlattenProcessGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in imageFlattenProcessGUI.M with the given input arguments.
%
%      imageFlattenProcessGUI('Property','Value',...) creates a new imageFlattenProcessGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imageFlattenProcessGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imageFlattenProcessGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imageFlattenProcessGUI

% Last Modified by GUIDE v2.5 03-Jul-2012 14:58:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imageFlattenProcessGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @imageFlattenProcessGUI_OutputFcn, ...
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


% --- Executes just before imageFlattenProcessGUI is made visible.
function imageFlattenProcessGUI_OpeningFcn(hObject, eventdata, handles, varargin)

processGUI_OpeningFcn(hObject, eventdata, handles, varargin{:},'initChannel',0);

% ---------------------- Channel Setup -------------------------
userData = get(handles.figure1, 'UserData');
funParams = userData.crtProc.funParams_;

% Set up available input channels
set(handles.listbox_availableChannels,'String',userData.MD.getChannelPaths(), ...
    'UserData',1:numel(userData.MD.channels_));

channelIndex = funParams.ChannelIndex;

% Find any parent process
userData.parentProc = userData.crtPackage.getParent(userData.procID);
if isempty(userData.crtPackage.processes_{userData.procID}) && ~isempty(userData.parentProc)
    % Check existence of all parent processes
    emptyParentProc = any(cellfun(@isempty,userData.crtPackage.processes_(userData.parentProc)));
    if ~emptyParentProc
        % Intersect channel index with channel index of parent processes
        parentChannelIndex = @(x) userData.crtPackage.processes_{x}.funParams_.ChannelIndex;
        for i = userData.parentProc
            channelIndex = intersect(channelIndex,parentChannelIndex(i));
        end
    end
   
end

if ~isempty(channelIndex)
    channelString = userData.MD.getChannelPaths(channelIndex);
else
    channelString = {};
end

set(handles.listbox_selectedChannels,'String',channelString,...
    'UserData',channelIndex);

set(handles.edit_GaussFilterSigma,'String',funParams.GaussFilterSigma);

% flattenMethods = userData.crtProc.getMethods();
% set(handles.popupmenu_flatteningMethods,'String',{flattenMethods(:).name},...
%     'Value',funParams.MethodIndx);

% Update user data and GUI data
handles.output = hObject;
set(hObject, 'UserData', userData);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = imageFlattenProcessGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% Notify the package GUI that the setting panel is closed
userData = get(handles.figure1, 'UserData');

set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);


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
% update_data(hObject,eventdata,handles);

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
% update_data(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Call back function of 'Apply' button
userData = get(handles.figure1, 'UserData');

% -------- Check user input --------
if isempty(get(handles.listbox_selectedChannels, 'String'))
   errordlg('Please select at least one input channel from ''Available Channels''.','Setting Error','modal') 
    return;
end
channelIndex = get (handles.listbox_selectedChannels, 'Userdata');
funParams.ChannelIndex = channelIndex;

gaussFilterSigma = str2double(get(handles.edit_GaussFilterSigma, 'String'));
if isnan(gaussFilterSigma) || gaussFilterSigma < 0
    errordlg(['Please provide a valid input for '''...
        get(handles.text_GaussFilterSigma,'String') '''.'],'Setting Error','modal');
    return;
end
funParams.GaussFilterSigma=gaussFilterSigma;

funParams.method_ind = 1; 

funParams.OutputDirectory  = [ userData.crtPackage.outputDirectory_, filesep 'ImageFlatten'];

for iChannel = channelIndex
ImageFlattenChannelOutputDir = [funParams.OutputDirectory,'/Channel',num2str(iChannel)];
    if (~exist(ImageFlattenChannelOutputDir,'dir'))
        mkdir(ImageFlattenChannelOutputDir);
    end
    
    userData.crtProc.setOutImagePath(iChannel,ImageFlattenChannelOutputDir)
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

% Set parameters and update main window
processGUI_ApplyFcn(hObject, eventdata, handles, funParams);


function update_data(hObject,eventdata, handles)

userData = get(handles.figure1, 'UserData');

% if strcmp(get(get(hObject,'Parent'),'Tag'),'uipanel_channels') ||...
%         strcmp(get(hObject,'Tag'),'listbox_thresholdValues')
%     % Check if changes have been at the list box level
%     linkedListBoxes = {'listbox_selectedChannels'};
%     checkLinkBox = strcmp(get(hObject,'Tag'),linkedListBoxes);
%     if any(checkLinkBox)
%         value = get(handles.(linkedListBoxes{checkLinkBox}),'Value');
%         set(handles.(linkedListBoxes{~checkLinkBox}),'Value',value);
%     else
%         value = get(handles.listbox_selectedChannels,'Value');
%     end  
% end


% % Retrieve the channex index
% props=get(handles.listbox_selectedChannels,{'UserData','Value'});
% if isempty(props{1}), return; end
% chanIndx = props{1}(props{2});
% imIndx = get(handles.slider_frameNumber,'Value');
% 
% % Load a new image in case the image number or channel has been changed
% if (chanIndx~=userData.chanIndx) ||  (imIndx~=userData.imIndx)
%     if ~isempty(userData.parentProc) && ~isempty(userData.crtPackage.processes_{userData.parentProc}) &&...
%             userData.crtPackage.processes_{userData.parentProc}.checkChannelOutput(chanIndx)
%         userData.imData=userData.crtPackage.processes_{userData.parentProc}.loadOutImage(chanIndx,imIndx);
%     else
%         userData.imData=userData.MD.channels_(chanIndx).loadImage(imIndx);
%     end
%     
%     userData.updateImage=1;
%     userData.chanIndx=chanIndx;
%     userData.imIndx=imIndx;
% else
%     userData.updateImage=0;
% end


% Save the data
set(handles.figure1, 'UserData', userData);
guidata(hObject,handles);
% 
% % Update graphics if applicable
% if get(handles.checkbox_preview,'Value')
% 
%     % Create figure if non-existing or closed
%     if ~ishandle(userData.previewFig)
%         userData.previewFig = figure('NumberTitle','off','Name','Thresholding preview',...
%             'Position',[50 50 userData.MD.imSize_(2) userData.MD.imSize_(1)]);
%         axes('Position',[0 0 1 1]);
%     else
%         figure(userData.previewFig);
%     end
%     
%     % Retrieve the handle of the figures image
%     imHandle = findobj(userData.previewFig,'Type','image');
%     if userData.updateImage || isempty(imHandle)
%         if isempty(imHandle)
%             imHandle=imagesc(userData.imData);
%             axis off;
%         else
%             set(imHandle,'CData',userData.imData);
%         end
%     end
%     
%     % Preview the tresholding output using the alphaData mapping    
%     flatten_preview_img=ones(size(userData.imData));
%  
%     
%     if get(handles.checkbox_auto,'Value')
%         methodIndx=get(handles.popupmenu_flatteningMethod,'Value');
%         
%         flatten_preview_img=image_flatten_single(currentImg, flatten_method_ind, Gaussian_sigma);
%     end
%     
%       gaussFilterSigma = str2double(get(handles.edit_GaussFilterSigma, 'String'));
%     if ~isnan(gaussFilterSigma) && gaussFilterSigma >0
%         flatten_preview_img = filterGauss2D(flatten_preview_img,gaussFilterSigma);    
%     end
%     
%     figure(imHandle); imagescc(flatten_preview_img);
%     
%     set(handles.figure1, 'UserData', userData);
%     guidata(hObject,handles);
% end

