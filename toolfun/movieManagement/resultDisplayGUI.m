function varargout = resultDisplayGUI(varargin)
% RESULTDISPLAYGUI M-file for resultDisplayGUI.fig
%      RESULTDISPLAYGUI, by itself, creates a new RESULTDISPLAYGUI or raises the existing
%      singleton*.
%
%      H = RESULTDISPLAYGUI returns the handle to a new RESULTDISPLAYGUI or the handle to
%      the existing singleton*.
%
%      RESULTDISPLAYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESULTDISPLAYGUI.M with the given input arguments.
%
%      RESULTDISPLAYGUI('Property','Value',...) creates a new RESULTDISPLAYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before resultDisplayGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to resultDisplayGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help resultDisplayGUI

% Last Modified by GUIDE v2.5 07-Oct-2010 17:08:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @resultDisplayGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @resultDisplayGUI_OutputFcn, ...
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


% --- Executes just before resultDisplayGUI is made visible.
function resultDisplayGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% resultDisplayGUI(process)
%
% User Data:
%   
%   userData.MD - movie data
%   userData.process - the process
%   userData.nFrames - total number of frames in this movie
%   userData.iFrame - current frame to display
%   userData.cMap - color map
%   userData.iChan - array of channel indexes who have valid result
%   userData.overlay - [], 'mask', 
%
userData = get(handles.figure1, 'UserData');
% Choose default command line output for recycleProcessGUI
handles.output = hObject;


if nargin > 3
    
    
    assert( isa(varargin{1}, 'Process'), 'User-defined: The first input must be Process object.')
    userData.process = varargin{1};
    
    % User Data
    userData.MD = userData.process.owner_;
    userData.nFrames = userData.MD.nFrames_;
    userData.iFrame = 1; % Default is to display the first frame
    userData.cMap = 'gray';
    userData.colStr = {'r','g','b'};
    userData.nChan = length(userData.MD.channels_);
    userData.iChan = find(userData.process.checkChannelOutput);
    userData.hText = [handles.text_1 handles.text_2 handles.text_3];
    
    userData.display = [ ];
    userData.overlay = [ ];
    userData.dispChan = zeros(1,3);
    userData.dispInd = 1;  % which layer to draw the next channel (rgb)
    
    userData.iFrame = 1;
    userData.hImage = [ ];
    userData.hOverlay_mask = cell(1, userData.nChan);
    
    if isempty(userData.iChan)
       return 
    end
    
    % GUI set-up
    set(handles.text_process, 'String', userData.process.name_)
    set(handles.text_movie, 'String', [userData.MD.movieDataPath_ userData.MD.movieDataFileName_])
    set(handles.edit_frame, 'String', num2str(userData.iFrame))
    set(handles.text_framenum, 'String', num2str(userData.nFrames))
    hold(handles.axes_1, 'on');
    
    % Create channels under channel menu
    handles.channel = arrayfun(@(x)uimenu(handles.menu_channel,'Label',['Channel ' num2str(x)],'Callback',@submenu_channel_Callback, 'Tag',['channel_' num2str(x)], 'UserData', x, 'Enable', 'off'), 1:userData.nChan );
	arrayfun(@(x)set(handles.channel(x), 'Enable', 'on'), userData.iChan)
    
    % Set display channels
    if length(userData.iChan) <= 3
        
        userData.dispChan(1:length(userData.iChan)) = userData.iChan;
        userData.dispInd = mod(length(userData.iChan), 3) + 1;

    else
        userData.dispChan(1:3) = userData.iChan(1:3);
        userData.dispInd = 1;

    end
    
    % Channel checked on
    arrayfun(@(x)set(handles.channel(x), 'Checked', 'on'), userData.dispChan(userData.dispChan > 0) )
    
    
% ---------------- Draw the image -------------------------------
    
    % Determine what type of process
    
    if isa(userData.process, 'ImageProcessingProcess') || isa(userData.process, 'MaskUtilityProcess')
        
        % Display Setting
        set(handles.menu_display_result, 'Enable', 'on', 'Checked', 'on')
        userData.display = 'menu_display_result';
        
        % Overlay Setting
        % Any?
        
    elseif isa(userData.process, 'SegmentationProcess')
        
        % Display Setting
        set(handles.menu_display_original, 'Checked', 'on')
        userData.display = 'menu_display_original';
        
        % Overlay Setting
        set(handles.menu_overlay_mask, 'Enable', 'on', 'Checked', 'on')
        userData.overlay = 'menu_overlay_mask';
        
    else
        error('User-defined: Not a valid process.')
    end
    
    set(handles.figure1,'UserData',userData)
    dispImage(userData.process, handles, userData.dispChan(userData.dispChan>0))
    userData = get(handles.figure1,'UserData');    
    
    % Display legend
    set(handles.figure1,'UserData',userData)
    
% ---------------- Draw the overlays -------------------------------

    if ~isempty(userData.overlay)
        switch userData.overlay

            case 'menu_overlay_mask'

                set(handles.figure1,'UserData',userData)
                dispOverlay_mask(userData.process, handles, userData.dispChan(userData.dispChan>0))
                userData = get(handles.figure1,'UserData');

        end
    end
    
end


set(handles.figure1,'UserData',userData)
dispLegend(handles)

guidata(hObject, handles)


function dispImage(process, handles, chanId, onoff)

if nargin < 4
    onoff = 1;
end

userData = get(handles.figure1, 'UserData');

% If input channel index is empty
if isempty(chanId)
%     set(handles.edit_notdisplay, 'Visible', 'on')
    return
end

% Check if all channel index are valid
if isempty(chanId) || any(chanId > userData.nChan)
   error('User-defined: channel index invalid') 
end

% Check if all channel index are all in userData.dispChan
if ~all( arrayfun(@(x)any(userData.dispChan == x), chanId) )
   error('User-defined: Not all input channel index are in the ') 
end


if onoff
% Show channels
    % Check if image has already existed
    if isempty(userData.hImage) || ~ishandle(userData.hImage)

        currImage = zeros([userData.MD.imSize_([2 1]) 3]);
    else

        currImage = get(userData.hImage, 'CData');
    end


    switch lower(userData.display)

        case 'menu_display_original'

            %Use the raw data as the image directories
            imDirs = userData.MD.getChannelPaths(chanId); 
            imNames = userData.MD.getImageFileNames(chanId); 

                for j = 1:length(chanId)
                    %load the image   
                    rgb = logical(userData.dispChan == chanId(j));

                    currImage(:,:,rgb) = mat2gray(imread([imDirs{j} filesep imNames{j}{userData.iFrame}]));  

                end
     

        case 'menu_display_result'  

                for j = 1:length(chanId)
                    %load the image    
                    rgb = logical(userData.dispChan == chanId(j));

                    currImage(:,:,rgb) = mat2gray(process.loadOutImage(chanId(j),userData.iFrame));   
                end

    end

    % Draw the image
    if isempty(userData.hImage) || ~ishandle(userData.hImage)

        set(handles.figure1, 'CurrentAxes', handles.axes_1)
        userData.hImage = imshow(currImage, []); 
    else
        set(userData.hImage, 'CData', currImage);
    end

else
% Hide channels
    if ~isempty(userData.hImage) && ishandle(userData.hImage)
        for j = 1:length(chanId)
            
            rgb = logical(userData.dispChan == chanId(j));
            currImage = get(userData.hImage, 'CData');
            currImage(:,:,rgb) = 0;
            set(userData.hImage, 'CData', currImage)            
        end
    end
    
    
end

set(handles.figure1, 'UserData', userData)


function dispOverlay_mask(segProcess, handles, chanId, onoff)
% Display mask overlay on image

if nargin < 4
   onoff = 1; 
end

userData = get(handles.figure1, 'UserData');

% If input channel index is empty
if isempty(chanId)
    disp('Warning: mask channel index is empty.')
    return
end

% Check if all channel index are valid
if any(chanId > userData.nChan)
   error('User-defined: channel index invalid') 
end

% Check if all channel index are all in userData.dispChan
if ~all( arrayfun(@(x)any(userData.dispChan == x), chanId) )
   error('User-defined: Not all input channel index are in the ') 
end

if onoff
% If turn mask on

    chanId = chanId( logical( segProcess.checkChannelOutput(chanId) ));
    
    if isempty(chanId)
       disp('Warning: No selected channel has mask output.') 
       return
    end

    for j = chanId
        
        rgb = logical(userData.dispChan == j);

%         if hasMasks(j) && userData.dispChan(j)>0

            if isempty(userData.hOverlay_mask{j})

                maskNames = segProcess.getOutMaskFileNames(j);

                %Load the mask
                currMask = imread([ segProcess.outMaskPaths_{j} filesep maskNames{1}{userData.iFrame}]);

                %Convert the mask into a boundary
                maskBounds = bwboundaries(currMask);

                %Plot the boundar(ies)                
                userData.hOverlay_mask{j} = cellfun(@(x)(plot(x(:,2),x(:,1),userData.colStr{rgb})),maskBounds); 
            else

                arrayfun( @(x)set(x, 'Visible', 'on', 'color', userData.colStr{rgb}), userData.hOverlay_mask{j} )
                
            end
%          end                                    
    end

else
% If turn mask off

    for j = chanId

        if  ~isempty(userData.hOverlay_mask{j})

            arrayfun( @(x)set(x, 'Visible', 'off'), userData.hOverlay_mask{j} )
        end

    end
    
end

set(handles.figure1, 'UserData', userData)            



function dispLegend(handles)
% GUI function: display the legend RGB baside image axes 

userData = get(handles.figure1, 'UserData');

%Draw the channel name
dispChan = userData.dispChan;

j = 0;
i = 0;
for k = 1 : length(dispChan)
    if dispChan(k) > 0
        
        j = j+1;
        set(userData.hText(j), 'Visible', 'on', 'ForegroundColor', userData.colStr{k}, ...
            'String', [upper(userData.colStr{k}) ': Channel ' num2str(dispChan(k))])
    else
        
        set(userData.hText(end-i), 'Visible', 'off')
        i = i +1;
    end

end
%Draw the frame number, or time in seconds if available
% if ~isempty(userData.MD.timeInterval_)
%     text(20,80,[num2str((userData.iFrame-1)*userData.MD.timeInterval_) ' / ' ...
%         num2str((userData.MD.nFrames_-1)*userData.MD.timeInterval_) ' s' ],'color',userData.colStr{1},'FontSize',12)
% else
%     text(20,80,[num2str(userData.iFrame) ' / ' num2str(userData.MD.nFrames_)],...
%         'color',userData.colStr{1},'FontSize',12)
% end


function dispNewFrame(handles)

userData = get(handles.figure1, 'UserData');

% Draw new frame
cla(handles.axes_1)
userData.hImage = [ ];
userData.hOverlay_mask = cell(1, userData.nChan);

set(handles.figure1, 'UserData', userData)
dispImage(userData.process, handles, userData.dispChan(userData.dispChan>0))

if ~isempty(userData.overlay)
    dispOverlay_mask(userData.process, handles, userData.dispChan(userData.dispChan>0))
end



% --- Outputs from this function are returned to the command line.
function varargout = resultDisplayGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if userData.iFrame == 1
    return
end

userData.iFrame = userData.iFrame - 1;
set(handles.edit_frame, 'String', num2str(userData.iFrame))
set(handles.figure1, 'UserData', userData)

if userData.iFrame == 1
    set(handles.pushbutton_left, 'Enable', 'off')
    set(handles.pushbutton_first, 'Enable', 'off')
end

if userData.iFrame == userData.nFrames - 1
   set(handles.pushbutton_right, 'Enable', 'on')
   set(handles.pushbutton_last, 'Enable', 'on')
end

% Display new frame
set(handles.figure1, 'UserData', userData)
dispNewFrame(handles)
    


% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if userData.iFrame == userData.nFrames
    return
end

userData.iFrame = userData.iFrame + 1;
set(handles.edit_frame, 'String', num2str(userData.iFrame))
set(handles.figure1, 'UserData', userData)

if userData.iFrame == userData.nFrames
    set(handles.pushbutton_right, 'Enable', 'off')
    set(handles.pushbutton_last, 'Enable', 'off')
end

if userData.iFrame == 2
   set(handles.pushbutton_left, 'Enable', 'on') 
   set(handles.pushbutton_first, 'Enable', 'on')
end
    
% Display new frame
set(handles.figure1, 'UserData', userData)
dispNewFrame(handles)



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_channel_Callback(hObject, eventdata, handles)
% hObject    handle to menu_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function submenu_channel_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
userData = get(handles.figure1, 'UserData');
channelid = get(hObject, 'UserData');

if strcmp(get(hObject, 'Checked'), 'on')
% Check off    
    
    % Get userData.dipChan
    zr = (userData.dispChan == channelid);
    id = find(zr);
    
    if ~any(zr)
        error('User-defined: no channel index found in userData.dispChan.') 
    elseif length(id)>1
        error('User-defined: more than one same channel index found in userData.dispChan.')
    end
    
    set(hObject, 'Checked', 'off')
    
    % Update CData
    if ~isfield(userData, 'hImage') || ~ishandle(userData.hImage)
       error('User-defined: no userData.hImage is found when updating CData.') 
    end
    
    CData = get(userData.hImage, 'CData');
    CData(:,:,id) = 0;
    set(userData.hImage, 'CData', CData)
    
    % Update Overlay
    if ~isempty(userData.overlay) && ~isempty(userData.hOverlay_mask{userData.dispChan(id)})
        arrayfun( @(x)set(x, 'Visible', 'off'), userData.hOverlay_mask{userData.dispChan(id)} )
    end
    
    % Update userData.dispChan
    userData.dispChan(id) = 0;   
    
    % Find next display indicator
    set(handles.figure1, 'UserData', userData)
    getNextDispInd(handles);
    userData = get(handles.figure1, 'UserData');
    
else
% Check on

    if any(userData.dispChan == channelid)
       error('User-defined: The channel index specified already exists in userData.dispChan.') 
    end
    
    set(hObject, 'Checked', 'on')
    
    % Uncheck old channel and hide its overlay
    if userData.dispChan(userData.dispInd) ~= 0
        
        set(handles.channel(userData.dispChan(userData.dispInd)), 'Checked', 'off')

        if ~isempty(userData.overlay) && ~isempty(userData.hOverlay_mask{userData.dispChan(userData.dispInd)})
            arrayfun( @(x)set(x, 'Visible', 'off'), userData.hOverlay_mask{userData.dispChan(userData.dispInd)} )
        end
    end
    
    % Update userData.dispChan
    userData.dispChan(userData.dispInd) = channelid;
    
    % Update CData
    set(handles.figure1, 'UserData', userData)
    dispImage(userData.process, handles, channelid, 1)
    userData = get(handles.figure1, 'UserData');    
    
    % Create new or show existing overlay
    if ~isempty(userData.overlay)

        set(handles.figure1, 'UserData', userData)
        dispOverlay_mask(userData.process, handles, channelid, 1)
        userData = get(handles.figure1, 'UserData');

    end
    
    % Find next display indicator
    set(handles.figure1, 'UserData', userData)
    getNextDispInd(handles)
    userData = get(handles.figure1, 'UserData');    
    
    
end

set(handles.figure1, 'UserData', userData)
dispLegend(handles)


function submenu_display_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if strcmp(get(hObject, 'Checked'), 'on')    
    return
end

% Check on new item
set(hObject, 'Checked', 'on')

switch get(hObject, 'tag')
    
    case 'menu_display_original'
        
        set(handles.menu_display_result, 'Checked', 'off')
        userData.display = 'menu_display_original';
          
    case 'menu_display_result'
        
        set(handles.menu_display_original, 'Checked', 'off')
        userData.display = 'menu_display_result';
        
    otherwise
        error('User-defined: internal error.')
        
end

set(handles.figure1,'UserData',userData)
dispImage(userData.process, handles, userData.dispChan(userData.dispChan>0))



function submenu_overlay_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if strcmp(get(hObject, 'Checked'), 'on')
% Check off    
    set(hObject, 'Checked', 'off')
    userData.overlay = [ ];
   
else
% Check on
    set(hObject, 'Checked', 'on')
    
    % Check off last image
    if ~isempty( userData.overlay )
        eval(['set(handles.' userData.overlay ', ''Checked'', ''off'')'])
    end    
    userData.overlay = get(hObject, 'tag');
    
end

set(handles.figure1,'UserData',userData)
dispOverlay_mask(userData.process, handles, userData.dispChan(userData.dispChan>0), strcmp(get(hObject, 'Checked'), 'on'))


function submenu_colormap_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if strcmp(get(hObject, 'Checked'), 'on')
% Check off    
    set(hObject, 'Checked', 'off')
    
else
% Check on
    set(hObject, 'Checked', 'on')
end

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)



function edit_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userData = get(handles.figure1, 'UserData');
value = str2double( get(handles.edit_frame, 'String') );

if isnan(value) || value <= 0 || floor(value) ~= ceil(value)
    errordlg('Please provide a valid input for frame index.');
    set(handles.edit_frame, 'String', num2str(userData.iFrame))
    return;
    
elseif value > userData.nFrames
    errordlg(['The frame index you entered is larger than the number of frames: ' num2str(userData.nFrames) '.']);
    set(handles.edit_frame, 'String', num2str(userData.iFrame))
    return;  
    
elseif value == userData.iFrame
    return

end

if userData.iFrame == userData.nFrames
    
    set(handles.pushbutton_right, 'Enable', 'on')
    set(handles.pushbutton_last, 'Enable', 'on')
    
elseif userData.iFrame == 1
    
    set(handles.pushbutton_left, 'Enable', 'on')
    set(handles.pushbutton_first, 'Enable', 'on')    
    
end

userData.iFrame = value;

if userData.iFrame == userData.nFrames
    
    set(handles.pushbutton_right, 'Enable', 'off')
    set(handles.pushbutton_last, 'Enable', 'off')
    
elseif userData.iFrame == 1
    
    set(handles.pushbutton_left, 'Enable', 'off')
    set(handles.pushbutton_first, 'Enable', 'off')    
    
end

% Display new frame
set(handles.figure1, 'UserData', userData)
dispNewFrame(handles)






% --- Executes during object creation, after setting all properties.
function edit_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_last.
function pushbutton_last_Callback(hObject, eventdata, handles)
userData = get(handles.figure1, 'UserData');

if userData.iFrame == userData.nFrames
    return
end

if userData.iFrame == 1
   set(handles.pushbutton_left, 'Enable', 'on')
   set(handles.pushbutton_first, 'Enable', 'on')
end

userData.iFrame = userData.nFrames;
set(handles.edit_frame, 'String', num2str(userData.iFrame))
set(handles.figure1, 'UserData', userData)

set(handles.pushbutton_right, 'Enable', 'off')
set(handles.pushbutton_last, 'Enable', 'off')

% Display new frame
set(handles.figure1, 'UserData', userData)
dispNewFrame(handles)



% --- Executes on button press in pushbutton_first.
function pushbutton_first_Callback(hObject, eventdata, handles)

userData = get(handles.figure1, 'UserData');

if userData.iFrame == 1
    return
end

if userData.iFrame == userData.nFrames
   set(handles.pushbutton_right, 'Enable', 'on')
   set(handles.pushbutton_last, 'Enable', 'on')
end

userData.iFrame = 1;
set(handles.edit_frame, 'String', num2str(userData.iFrame))
set(handles.figure1, 'UserData', userData)

set(handles.pushbutton_left, 'Enable', 'off')
set(handles.pushbutton_first, 'Enable', 'off')

% Display new frame
set(handles.figure1, 'UserData', userData)
dispNewFrame(handles)

% --------------------------------------------------------------------
function menu_display_Callback(hObject, eventdata, handles)
% hObject    handle to menu_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_overlay_Callback(hObject, eventdata, handles)
% hObject    handle to menu_overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Channel bar clicked!')

% --------------------------------------------------------------------
function menu_colormap_Callback(hObject, eventdata, handles)
% hObject    handle to menu_colormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Channel bar clicked!')

% --------------------------------------------------------------------
function menu_overlay_mask_Callback(hObject, eventdata, handles)
% hObject    handle to menu_overlay_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_overlay_windows_Callback(hObject, eventdata, handles)
% hObject    handle to menu_overlay_windows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_overlay_protrusion_Callback(hObject, eventdata, handles)
% hObject    handle to menu_overlay_protrusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_display_original_Callback(hObject, eventdata, handles)
% hObject    handle to menu_display_original (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_display_result_Callback(hObject, eventdata, handles)
% hObject    handle to menu_display_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function getNextDispInd(handles)

userData = get(handles.figure1, 'UserData');

if all(userData.dispChan)
    
    userData.dispInd = mod(userData.dispInd, 3) +1;
    
else
    
    zr = find(userData.dispChan == 0);
    userData.dispInd = zr(1);
end
set(handles.figure1, 'UserData', userData)


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)

if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
if strcmp(eventdata.Key, 'leftarrow')
    pushbutton_left_Callback(handles.pushbutton_left, [], handles);
end
if strcmp(eventdata.Key, 'rightarrow')
    pushbutton_right_Callback(handles.pushbutton_right, [], handles);
end
